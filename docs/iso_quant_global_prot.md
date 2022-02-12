# Isobaric Quantification: Proteomics {#iso-global}

This pipeline shows how to process global proteomics TMT data with `PlexedPiper`. We will use data package 3442, which is "PlexedPiperTestData global". In addition to PlexedPiper, we will also need MSnID (the basis for PlexedPiper) and PNNL.DMS.utils to interface with PNNL's DMS.




```r
## Install missing packages
if (!require("remotes", quietly = T)) install.packages("remotes")
git_packages <- c("MSnID@pnnl-master", "PlexedPiper", "PNNL.DMS.utils")
for (pkg_i in git_packages) {
  if (!require(sub("@.*", "", pkg_i), quietly = T, character.only = T))
    remotes::install_github(file.path("PNNL-Comp-Mass-Spec", pkg_i))
}
## ------------------------

library(MSnID)
library(PlexedPiper)
library(PNNL.DMS.utils)
```



The pipeline can be broken up into four major parts: prepare MS/MS identifications, prepare reporter ion intensities, create study design tables, and create a quantitative cross-tab. There is another step that is required for statistical testing, which is to create an MSnSet.

## Prepare MS/MS Identifications

### Read MS-GF+ Data

The first step in the preparation of the MS/MS identifications is to fetch the data. This can either be obtained from PNNL's DMS or from a local folder.  If working with the DMS, use `PNNL.DMS.utils::read_msgf_data_from_DMS`; otherwise, use `PlexedPiper::read_msgf_data`.


```r
## Get MS-GF+ results from local folder - not run
# Get file path
path_to_MSGF_results <- "path_to_msgf_results"
# Read MS-GF+ data from path
msnid <- read_msgf_data(path_to_MSGF_results)
```


```r
## Get MS-GF+ results from DMS
data_package_num <- 3442 # global proteomics
msnid <- read_msgf_data_from_DMS(data_package_num) # global proteomics
```

Normally, this would display a progress bar in the console as the data is being fetched. However, the output was suppressed to save space. We can view a summary of the MSnID object with the `show()` function.


```r
show(msnid)
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  48 
## #PSMs: 1156754 at 31 % FDR
## #peptides: 511617 at 61 % FDR
## #accessions: 128378 at 98 % FDR
```

This summary tells us that `msnid` consists of 4 spectrum files (datasets), and contains a total of 1,156,754 peptide-spectrum-matches (PSMs), 511,617 total peptides, and 128,378 total accessions (proteins). The reported FDR is the empirical **false-discovery rate**, which is calculated as the ratio of the number of unique decoy to unique non-decoy PSMs, peptides, or accessions.

### Correct Isotope Selection Error 

Carbon has two stable isotopes: $^{12}\text{C}$ and $^{13}\text{C}$, with natural abundances of 98.93% and 1.07%, respectively [@berglund_isotopic_2011]. That is, we expect that about 1 out of every 100 carbon atoms is naturally going to be a $^{13}\text{C}$, while the rest are $^{12}\text{C}$. In larger peptides with many carbon atoms, it is more likely that at least one atom will be a $^{13}\text{C}$ than all atoms will be $^{12}\text{C}$. In cases such as these, a non-monoisotopic ion will be selected by the instrument for fragmentation.

<div class="figure" style="text-align: center">
<img src="images/MS1_non_monoisotopic.PNG" alt="MS1 spectra with peak at non-monoisotopic precursor ion." width="75%" />
<p class="caption">(\#fig:MS1-peak)MS1 spectra with peak at non-monoisotopic precursor ion.</p>
</div>

In Figure \@ref(fig:MS1-peak), the monoisotopic ion (m/z of 1427.29) is not the most abundant, so it is not selected as the precursor. Instead, the ion with a $^{13}\text{C}$ in place of a $^{12}\text{C}$ is selected for fragmentation. We calculate the mass difference between these two ions as the difference between the mass-to-charge ratios multiplied by the ion charge. In this case, the mass difference is 1 Dalton, or about the difference between $^{13}\text{C}$ and $^{12}\text{C}$. (More accurately, the difference between these isotopes is 1.0033548378 Da.) While MS-GF+ is still capable of correctly identifying these peptides, the downstream calculations of mass measurement error need to be fixed because they are used for filtering later on (Section \@ref(global-peptide-filter)). The `correct_peak_selection` function corrects these mass measurement errors, and Figure \@ref(fig:mass-to-charge-diff) shows the distribution of the absolute mass measurement errors (in PPM) before and after correction.




```r
# Correct for isotope selection error
msnid <- correct_peak_selection(msnid)
```

<div class="figure" style="text-align: center">
<img src="iso_quant_global_prot_files/figure-html/mass-to-charge-diff-1.png" alt="Histogram of mass measurement errors before and after correction." width="75%" />
<p class="caption">(\#fig:mass-to-charge-diff)Histogram of mass measurement errors before and after correction.</p>
</div>


### Remove Contaminants 

Now, we will remove contaminants such as the trypsin that was used for protein digestion. We can see which contaminants will be removed with `accessions(msnid)[grepl("Contaminant", accessions(msnid))]`. To remove contaminants, we use `apply_filter` with an appropriate string that tells the function what rows to keep. In this case, we keep rows where the accession does not contain "Contaminant". We will use `show` to see how the counts change.


```r
# Remove contaminants
msnid <- apply_filter(msnid, "!grepl('Contaminant', accession)")
show(msnid)
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  48 
## #PSMs: 1155442 at 31 % FDR
## #peptides: 511196 at 61 % FDR
## #accessions: 128353 at 98 % FDR
```

We can see that the number of PSMs decreased by about 1300, peptides by ~400, and proteins by 25.

### MS/MS ID Filter: Peptide Level {#global-peptide-filter}

The next step is to filter the MS/MS identifications such that the empirical peptide-level FDR is less than some threshold and the number of identifications is maximized. We will use the $-log_{10}$ of the `PepQValue` column as one of our filtering criteria and assign it to a new column in `psms(msnid)` called `msmsScore`. The `PepQValue` column is the MS-GF+ Spectrum E-value, which reflects how well the theoretical and experimental fragmentation spectra match; therefore, high values of `msmsScore` indicate a good match (see Figure \@ref(fig:plot-msmsScore)).

<div class="figure" style="text-align: center">
<img src="iso_quant_global_prot_files/figure-html/plot-msmsScore-1.png" alt="Density plot of msmsScore." width="75%" />
<p class="caption">(\#fig:plot-msmsScore)Density plot of msmsScore.</p>
</div>
</br>

The other filtering criteria is the absolute deviation of the mass measurement error of the precursor ions in parts-per-million (ppm), which is assigned to the `absParentMassErrorPPM` column in `psms(msnid)` (see Figure \@ref(fig:plot-mass-error)).

<div class="figure" style="text-align: center">
<img src="iso_quant_global_prot_files/figure-html/plot-mass-error-1.png" alt="Density plot of absParentMassErrorPPM." width="75%" />
<p class="caption">(\#fig:plot-mass-error)Density plot of absParentMassErrorPPM.</p>
</div>

</br>

These new columns `msmsScore` and `absParentMassErrorPPM` are generated automatically by `filter_msgf_data`, so we don't need to worry about creating them ourselves.


```r
# 1% FDR filter at the peptide level
msnid <- filter_msgf_data(msnid, level = "peptide", fdr.max = 0.01)
show(msnid)
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  48 
## #PSMs: 471295 at 0.46 % FDR
## #peptides: 96663 at 1 % FDR
## #accessions: 27098 at 9.1 % FDR
```

We can see that filtering drastically reduces the number of PSMs, and the empirical peptide-level FDR is now 1%. However, notice that the empirical protein-level FDR is still fairly high.

### MS/MS ID Filter: Protein Level 

Now, we need to filter proteins so that the FDR is at most 1%. For each protein, we divide the number of associated peptides by its length and multiply this value by 1000. This new `peptides_per_1000aa` column is used as the filter criteria (Figure \@ref(fig:plot-num-pep)).

We will need the lengths of each protein, which can be obtained from the FASTA (pronounced FAST-AYE) file that contains the protein sequences used in the database search. The first three entries of the FASTA file are shown in Figure \@ref(fig:fasta-ex).

<div class="figure" style="text-align: center">
<img src="images/FASTA_example_MoTrPAC.PNG" alt="First three entries of the FASTA file." width="75%" />
<p class="caption">(\#fig:fasta-ex)First three entries of the FASTA file.</p>
</div>

</br>

The path to the FASTA file can be specified as a local file path or it can be obtained with `PNNL.DMS.utils::path_to_FASTA_used_by_DMS`. We will use the latter method.


```r
## Get path to FASTA file from local folder - not run
path_to_FASTA <- "some_folder/name_of_fasta_file.fasta"
```


```r
## Get path to FASTA file from DMS
path_to_FASTA <- path_to_FASTA_used_by_DMS(data_package_num)
```


```r
# Compute number of peptides per 1000 amino acids
msnid <- compute_num_peptides_per_1000aa(msnid, path_to_FASTA)
```

<div class="figure" style="text-align: center">
<img src="iso_quant_global_prot_files/figure-html/plot-num-pep-1.png" alt="Density plot of peptides_per_1000aa. The plot area has been zoomed in." width="75%" />
<p class="caption">(\#fig:plot-num-pep)Density plot of peptides_per_1000aa. The plot area has been zoomed in.</p>
</div>
</br>

Now, we filter the proteins to 1% FDR.


```r
# 1% FDR filter at the protein level
msnid <- filter_msgf_data(msnid, level = "accession", fdr.max = 0.01)
show(msnid)
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  48 
## #PSMs: 464788 at 0.17 % FDR
## #peptides: 92176 at 0.34 % FDR
## #accessions: 15629 at 0.98 % FDR
```


### Inference of Parsimonious Protein Set 

The situation when a certain peptide sequence matches multiple proteins adds complication to the downstream quantitative analysis, as it is not clear which protein this peptide is originating from. There are common ways for dealing with this. One is to simply retain uniquely matching peptides and discard shared peptides (`unique_only = TRUE`). Alternatively, assign the shared peptides to the proteins with the larger number of uniquely mapping peptides (`unique_only = FALSE`). If there is a choice between multiple proteins with equal numbers of uniquely mapping peptides, the shared peptides are assigned to the first protein according to alphanumeric order (Figure \@ref(fig:parsimony)). 

<!---
This step could be done prior to filtering at the accession level, but if peptides are assigned to a low-confidence protein, and that protein is removed during filtering, those peptides will be lost. Instead, it is better to filter to the set of confidently-identified proteins and then determine the parsimonious set.
--->

<div class="figure" style="text-align: center">
<img src="images/parsimonious-protein-set-inference.PNG" alt="Visual explanation of the inference of the parsimonious protein set." width="75%" />
<p class="caption">(\#fig:parsimony)Visual explanation of the inference of the parsimonious protein set.</p>
</div>

</br>


```r
# Inference of parsimonious protein set
msnid <- infer_parsimonious_accessions(msnid, unique_only = FALSE)
show(msnid)
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  48 
## #PSMs: 451568 at 0.16 % FDR
## #peptides: 90639 at 0.28 % FDR
## #accessions: 5247 at 1.1 % FDR
```

Notice that the protein-level FDR increased slightly above the 1% threshold. In this case, the difference isn't significant, so we can ignore it. 

**Note:** If the peptide or accession-level FDR increases significantly above 1% after inference of the parsimonious protein set, consider lowering the FDR cutoff (for example, to 0.9%) and redoing the previous processing steps. That is, start with the MSnID prior to any filtering and redo the FDR filtering steps.

### Remove Decoy PSMs

The final step in preparing the MS/MS identifications is to remove the decoy PSMs, as they were only needed for the FDR filters. We use the `apply_filter` function again and only keep entries where `isDecoy` is `FALSE`.


```r
# Remove Decoy PSMs
msnid <- apply_filter(msnid, "!isDecoy")
show(msnid)
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  48 
## #PSMs: 450857 at 0 % FDR
## #peptides: 90382 at 0 % FDR
## #accessions: 5191 at 0 % FDR
```

After processing, we are left with 450,928 PSMs, 90,411 peptides, and 5,201 proteins. Table \@ref(tab:global-msnid-table) shows the first 6 rows of the processed MS-GF+ output.

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:90%; "><table class="table table-hover table-condensed" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">(\#tab:global-msnid-table)<left>First 6 rows of the processed MS-GF+ results.</left>
</caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Dataset </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ResultID </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Scan </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> FragMethod </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SpecIndex </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Charge </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PrecursorMZ </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> DelM </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> DelM_PPM </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> MH </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> peptide </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Protein </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> NTT </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> DeNovoScore </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> MSGFScore </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> MSGFDB_SpecEValue </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Rank_MSGFDB_SpecEValue </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> EValue </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> QValue </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PepQValue </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> IsotopeError </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> accession </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> calculatedMassToCharge </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> chargeState </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> experimentalMassToCharge </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> isDecoy </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> spectrumFile </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> spectrumID </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> pepSeq </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> msmsScore </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> absParentMassErrorPPM </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> peptides_per_1000aa </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_07_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 1862 </td>
   <td style="text-align:right;"> 27707 </td>
   <td style="text-align:left;"> HCD </td>
   <td style="text-align:right;"> 324 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 928.541 </td>
   <td style="text-align:right;"> -0.001 </td>
   <td style="text-align:right;"> -0.526 </td>
   <td style="text-align:right;"> 1856.075 </td>
   <td style="text-align:left;"> R.AAAAAAAAAAAAAAGAAGK.E </td>
   <td style="text-align:left;"> NP_113986.1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 285 </td>
   <td style="text-align:right;"> 282 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> NP_113986.1 </td>
   <td style="text-align:right;"> 928.541 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 928.541 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_07_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 27707 </td>
   <td style="text-align:left;"> AAAAAAAAAAAAAAGAAGK </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> 0.589 </td>
   <td style="text-align:right;"> 24.938 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_07_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 4192 </td>
   <td style="text-align:right;"> 27684 </td>
   <td style="text-align:left;"> HCD </td>
   <td style="text-align:right;"> 906 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 619.363 </td>
   <td style="text-align:right;"> -0.002 </td>
   <td style="text-align:right;"> -0.887 </td>
   <td style="text-align:right;"> 1856.075 </td>
   <td style="text-align:left;"> R.AAAAAAAAAAAAAAGAAGK.E </td>
   <td style="text-align:left;"> NP_113986.1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 156 </td>
   <td style="text-align:right;"> 144 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> NP_113986.1 </td>
   <td style="text-align:right;"> 619.363 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 619.363 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_07_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 27684 </td>
   <td style="text-align:left;"> AAAAAAAAAAAAAAGAAGK </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> 0.991 </td>
   <td style="text-align:right;"> 24.938 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_06_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 26263 </td>
   <td style="text-align:right;"> 27336 </td>
   <td style="text-align:left;"> HCD </td>
   <td style="text-align:right;"> 5187 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 619.363 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.197 </td>
   <td style="text-align:right;"> 1856.075 </td>
   <td style="text-align:left;"> R.AAAAAAAAAAAAAAGAAGK.E </td>
   <td style="text-align:left;"> NP_113986.1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 118 </td>
   <td style="text-align:right;"> 85 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> NP_113986.1 </td>
   <td style="text-align:right;"> 619.363 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 619.363 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_06_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 27336 </td>
   <td style="text-align:left;"> AAAAAAAAAAAAAAGAAGK </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> 0.091 </td>
   <td style="text-align:right;"> 24.938 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_07_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 1471 </td>
   <td style="text-align:right;"> 27096 </td>
   <td style="text-align:left;"> HCD </td>
   <td style="text-align:right;"> 415 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 619.363 </td>
   <td style="text-align:right;"> -0.001 </td>
   <td style="text-align:right;"> -0.591 </td>
   <td style="text-align:right;"> 1856.075 </td>
   <td style="text-align:left;"> R.AAAAAAAAAAAAAAGAAGK.E </td>
   <td style="text-align:left;"> NP_113986.1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 157 </td>
   <td style="text-align:right;"> 156 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> NP_113986.1 </td>
   <td style="text-align:right;"> 619.363 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 619.363 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_07_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 27096 </td>
   <td style="text-align:left;"> AAAAAAAAAAAAAAGAAGK </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> 0.684 </td>
   <td style="text-align:right;"> 24.938 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_05_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 28664 </td>
   <td style="text-align:right;"> 10441 </td>
   <td style="text-align:left;"> HCD </td>
   <td style="text-align:right;"> 4849 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 586.832 </td>
   <td style="text-align:right;"> -0.001 </td>
   <td style="text-align:right;"> -0.728 </td>
   <td style="text-align:right;"> 1172.659 </td>
   <td style="text-align:left;"> R.AAAAADLANR.S </td>
   <td style="text-align:left;"> NP_001007804.1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 124 </td>
   <td style="text-align:right;"> 124 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.002 </td>
   <td style="text-align:right;"> 0.003 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> NP_001007804.1 </td>
   <td style="text-align:right;"> 586.833 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 586.832 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_05_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 10441 </td>
   <td style="text-align:left;"> AAAAADLANR </td>
   <td style="text-align:right;"> 2.480 </td>
   <td style="text-align:right;"> 0.746 </td>
   <td style="text-align:right;"> 34.755 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_24_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 41775 </td>
   <td style="text-align:right;"> 8033 </td>
   <td style="text-align:left;"> HCD </td>
   <td style="text-align:right;"> 7889 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 831.447 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 1661.886 </td>
   <td style="text-align:left;"> G.AAAAAEAESGGGGGK.K </td>
   <td style="text-align:left;"> NP_001128630.1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 176 </td>
   <td style="text-align:right;"> 76 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.001 </td>
   <td style="text-align:right;"> 0.003 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> NP_001128630.1 </td>
   <td style="text-align:right;"> 831.447 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 831.447 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_24_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 8033 </td>
   <td style="text-align:left;"> AAAAAEAESGGGGGK </td>
   <td style="text-align:right;"> 2.583 </td>
   <td style="text-align:right;"> 0.106 </td>
   <td style="text-align:right;"> 579.815 </td>
  </tr>
</tbody>
</table></div>

</br>


## Prepare Reporter Ion Intensities

### Read MASIC Output 

MASIC is a tool for extracting ion intensities. With proper parameter settings, it can be used for extracting TMT (or iTRAQ) reporter ion intensities. In addition, it reports a number of other helpful metrics. Notably, the interference score at the precursor ion level and the signal-to-noise ratio (S/N) at the reporter ion level (computed by Thermo software). The interference score reflects the proportion of the ion population that was isolated for fragmentation that is due to the targeted ion. In other words, `1 - InterferenceScore` is due to co-isolated species that have similar elution time and precursor ion m/z. The first step in the preparation of the reporter ion intensity data is to read the MASIC results. By default, the interference score is not included, so we need to set that argument to `TRUE` in order to filter the results.

Similar to the MS-GF+ results, we can read the MASIC results from a local folder with `PlexedPiper::read_masic_data` or from PNNL's DMS with `PNNL.DMS.utils::read_masic_data_from_DMS`.


```r
## Get MASIC results from local folder - not run
# Get file path
path_to_MASIC_results <- "path_to_folder_containing_necessary_files"
# Read MASIC results from path
masic_data <- read_masic_data(path_to_MASIC_results, 
                              interference_score = TRUE)
```


```r
## Get MASIC results from DMS
masic_data <- read_masic_data_from_DMS(data_package_num, 
                                       interference_score = TRUE)
```

Normally, this would display progress bars in the console as the data is being fetched. However, the output was suppressed to save space.

Table \@ref(tab:global-masic-unfiltered) shows the first 6 rows of the unfiltered `masic_data`.

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:90%; "><table class="table table-hover table-condensed" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">(\#tab:global-masic-unfiltered)<left>First 6 rows of the unfiltered MASIC data.</left>
</caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Dataset </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ScanNumber </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_126.128 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_127.125 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_127.131 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_128.128 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_128.134 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_129.131 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_129.138 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_130.135 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_130.141 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_131.138 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_126.128_SignalToNoise </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_127.125_SignalToNoise </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_127.131_SignalToNoise </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_128.128_SignalToNoise </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_128.134_SignalToNoise </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_129.131_SignalToNoise </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_129.138_SignalToNoise </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_130.135_SignalToNoise </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_130.141_SignalToNoise </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_131.138_SignalToNoise </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> InterferenceScore </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 70562.39 </td>
   <td style="text-align:right;"> 24864.62 </td>
   <td style="text-align:right;"> 17165.80 </td>
   <td style="text-align:right;"> 35625.00 </td>
   <td style="text-align:right;"> 92236.87 </td>
   <td style="text-align:right;"> 9640.23 </td>
   <td style="text-align:right;"> 8578.05 </td>
   <td style="text-align:right;"> 6996.69 </td>
   <td style="text-align:right;"> 11833.07 </td>
   <td style="text-align:right;"> 32281.34 </td>
   <td style="text-align:right;"> 71.47 </td>
   <td style="text-align:right;"> 25.17 </td>
   <td style="text-align:right;"> 17.38 </td>
   <td style="text-align:right;"> 36.04 </td>
   <td style="text-align:right;"> 93.32 </td>
   <td style="text-align:right;"> 9.75 </td>
   <td style="text-align:right;"> 8.67 </td>
   <td style="text-align:right;"> 7.07 </td>
   <td style="text-align:right;"> 11.96 </td>
   <td style="text-align:right;"> 32.71 </td>
   <td style="text-align:right;"> 0.996 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 23706.89 </td>
   <td style="text-align:right;"> 13559.32 </td>
   <td style="text-align:right;"> 5856.83 </td>
   <td style="text-align:right;"> 16322.71 </td>
   <td style="text-align:right;"> 34294.90 </td>
   <td style="text-align:right;"> 4853.11 </td>
   <td style="text-align:right;"> 7938.24 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 1465.03 </td>
   <td style="text-align:right;"> 18182.27 </td>
   <td style="text-align:right;"> 26.12 </td>
   <td style="text-align:right;"> 14.94 </td>
   <td style="text-align:right;"> 6.45 </td>
   <td style="text-align:right;"> 17.97 </td>
   <td style="text-align:right;"> 37.77 </td>
   <td style="text-align:right;"> 5.34 </td>
   <td style="text-align:right;"> 8.74 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 1.61 </td>
   <td style="text-align:right;"> 19.93 </td>
   <td style="text-align:right;"> 0.993 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 12459.86 </td>
   <td style="text-align:right;"> 11785.91 </td>
   <td style="text-align:right;"> 10932.51 </td>
   <td style="text-align:right;"> 10653.32 </td>
   <td style="text-align:right;"> 12328.62 </td>
   <td style="text-align:right;"> 5959.86 </td>
   <td style="text-align:right;"> 9905.82 </td>
   <td style="text-align:right;"> 8387.04 </td>
   <td style="text-align:right;"> 11166.70 </td>
   <td style="text-align:right;"> 14053.40 </td>
   <td style="text-align:right;"> 12.40 </td>
   <td style="text-align:right;"> 11.75 </td>
   <td style="text-align:right;"> 10.90 </td>
   <td style="text-align:right;"> 10.64 </td>
   <td style="text-align:right;"> 12.31 </td>
   <td style="text-align:right;"> 5.96 </td>
   <td style="text-align:right;"> 9.91 </td>
   <td style="text-align:right;"> 8.40 </td>
   <td style="text-align:right;"> 11.18 </td>
   <td style="text-align:right;"> 14.13 </td>
   <td style="text-align:right;"> 1.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 1.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 10998.67 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 21077.05 </td>
   <td style="text-align:right;"> 2725.50 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 6800.70 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 9.19 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 17.57 </td>
   <td style="text-align:right;"> 2.27 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 5.66 </td>
   <td style="text-align:right;"> 1.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 6166.82 </td>
   <td style="text-align:right;"> 1371.27 </td>
   <td style="text-align:right;"> 2418.35 </td>
   <td style="text-align:right;"> 8087.76 </td>
   <td style="text-align:right;"> 5485.35 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 1543.48 </td>
   <td style="text-align:right;"> 1943.96 </td>
   <td style="text-align:right;"> 7436.60 </td>
   <td style="text-align:right;"> 6.92 </td>
   <td style="text-align:right;"> 1.54 </td>
   <td style="text-align:right;"> 2.71 </td>
   <td style="text-align:right;"> 9.04 </td>
   <td style="text-align:right;"> 6.13 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 1.72 </td>
   <td style="text-align:right;"> 2.16 </td>
   <td style="text-align:right;"> 8.26 </td>
   <td style="text-align:right;"> 0.969 </td>
  </tr>
</tbody>
</table></div>

</br>


### Filter MASIC Data 

The only other step in reporter ion intensity data preparation is to filter the results. Currently, we recommend keeping entries where at least 50% of the ion population is due to the targeted ion (interference score $\geq$ 0.5) and not filtering by S/N.


```r
# Filter MASIC data
masic_data <- filter_masic_data(masic_data, 
                                interference_score_threshold = 0.5,
                                s2n_threshold = 0)
```

Table \@ref(tab:global-masic-unfiltered) shows the first 6 rows of the filtered `masic_data`.

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:90%; "><table class="table table-hover table-condensed" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">(\#tab:global-masic-filtered)<left>First 6 rows of the filtered MASIC data.</left>
</caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Dataset </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ScanNumber </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_126.128 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_127.125 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_127.131 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_128.128 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_128.134 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_129.131 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_129.138 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_130.135 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_130.141 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_131.138 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 70562.39 </td>
   <td style="text-align:right;"> 24864.62 </td>
   <td style="text-align:right;"> 17165.80 </td>
   <td style="text-align:right;"> 35625.00 </td>
   <td style="text-align:right;"> 92236.87 </td>
   <td style="text-align:right;"> 9640.23 </td>
   <td style="text-align:right;"> 8578.05 </td>
   <td style="text-align:right;"> 6996.69 </td>
   <td style="text-align:right;"> 11833.07 </td>
   <td style="text-align:right;"> 32281.34 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 23706.89 </td>
   <td style="text-align:right;"> 13559.32 </td>
   <td style="text-align:right;"> 5856.83 </td>
   <td style="text-align:right;"> 16322.71 </td>
   <td style="text-align:right;"> 34294.90 </td>
   <td style="text-align:right;"> 4853.11 </td>
   <td style="text-align:right;"> 7938.24 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 1465.03 </td>
   <td style="text-align:right;"> 18182.27 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 12459.86 </td>
   <td style="text-align:right;"> 11785.91 </td>
   <td style="text-align:right;"> 10932.51 </td>
   <td style="text-align:right;"> 10653.32 </td>
   <td style="text-align:right;"> 12328.62 </td>
   <td style="text-align:right;"> 5959.86 </td>
   <td style="text-align:right;"> 9905.82 </td>
   <td style="text-align:right;"> 8387.04 </td>
   <td style="text-align:right;"> 11166.70 </td>
   <td style="text-align:right;"> 14053.40 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 10998.67 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 21077.05 </td>
   <td style="text-align:right;"> 2725.50 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 6800.70 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 6166.82 </td>
   <td style="text-align:right;"> 1371.27 </td>
   <td style="text-align:right;"> 2418.35 </td>
   <td style="text-align:right;"> 8087.76 </td>
   <td style="text-align:right;"> 5485.35 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 1543.48 </td>
   <td style="text-align:right;"> 1943.96 </td>
   <td style="text-align:right;"> 7436.60 </td>
  </tr>
</tbody>
</table></div>

Lastly, we will save the processed MSnID and MASIC data to an .RData file with compression. This is useful in case we want to create different cross-tabs with new study design tables later on.


```r
# Save processed MSnID and MASIC data
save(msnid, masic_data, file = "data/3442_processed_msnid_and_masic.RData",
     compress = TRUE)
```


## Create Study Design Tables {#fetch-study-design-tables}

To convert from PSMs and reporter ion intensities to meaningful quantitative data, it is necessary to specify the study design. The entire study design is captured by three tables - fractions, samples, references. With newly processed data, these typically do not exist, and must be created. The next sections show how to create these tables in R. 

**NOTE:** simple study designs can be created in Excel and read in with `readxl::read_excel`, though R is the better choice when dealing with many samples.


### Fractions 

The fractions table consists of two columns: `Dataset` and `PlexID`. The `Dataset` column contains all of the unique datasets that are common to `msnid` and `masic_data`. Sometimes, entire datasets may be removed during the FDR filtering steps, so that is why we use the unique intersection of datasets. The `PlexID` column contains the plex ID associated with each dataset, and is typically a letter followed by a number ("S1", "S2", etc.). A plex is a set of samples that are processed together (under the same conditions). We can extract the plex ID from the datasets. In this case, the plex ID always comes after "_W_", so we can use a regular expression (use `help(topic = regex, package = base)` to learn more). The regular expression below says to capture an "S" followed by a single digit that appears after "\_W\_" and before an underscore. The plex ID is always included in the dataset names, but the format of the names will be different.


```r
# Create fractions table
datasets <- unique(intersect(msnid$Dataset, masic_data$Dataset))
fractions <- data.frame(Dataset = datasets) %>% 
  mutate(PlexID = gsub(".*_W_(S\\d{1})_.*", "\\1", Dataset))
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:20em; "><table class="table table-hover table-condensed" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">(\#tab:fractions-table)<left>Fractions</left>
</caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Dataset </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> PlexID </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_07_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_06_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_07_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_05_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_24_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_08_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_16_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_16_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_23_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_15_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_09_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_17_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_08_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_23_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_22_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_24_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_05_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_04_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_04_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_02_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_13_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_12_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_03_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_03_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_18_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_18_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_20_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_06_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_21_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_10_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_19_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_15_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_17_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_12_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_14_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_13_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_14_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_09_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_11_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_22_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_10_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_11_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_19_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_20_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_21_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_02_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S2 </td>
  </tr>
</tbody>
</table></div>


### Samples 

The samples table contains columns `PlexID`, `QuantBlock`, `ReporterName`, `ReporterAlias`, and `MeasurementName`. 

- `PlexID` must be the same as the `PlexID` in the fractions table.
- `ReporterName` is the reporter ion name ("126", "127N", "127C", etc.). 
- `ReporterAlias` is used for defining the reference channel(s).
- `MeasurementName` determines the column names for the final cross-tab, and must be unique and begin with a letter. If any values of `ReporterAlias` are "ref", the corresponding `MeasurementName` should be `NA`. `NA` measurement names will not appear as columns in the final cross-tab.
- `QuantBlock` defines the sub-plex. In a typical TMT experiment, `QuantBlock` is always 1. In case of 5 pairwise comparisons within TMT10, there will be 5 QuantBlocks (1-5) with a reference for each `QuantBlock`.

For this experiment, TMT10 was used as the basis for two plexes, and channel 131 is the reference, so we set `ReporterAlias` to "ref" and `MeasurementName` to `NA` when `ReporterName` is `"131"`. This will divide the intensities of each channel by their associated reference and make the reference channel absent from the quantitative cross-tab. In cases where reporter ion intensities are not normalized by a reference channel (reference = 1) or they are normalized by the average of select channels, do not set any `ReporterAlias` to "ref" or `MeasurementName` to `NA`.


```r
# Create samples table
samples <- reporter_converter$tmt10 %>% 
  dplyr::select(ReporterName) %>% # only keep ReporterName column
  dplyr::slice(rep(1:n(), times = 2)) %>% # Copy TMT10 table twice (2 plexes)
  # Create PlexID and QuantBlock columns. 
  # Plex S1 goes with first 10 rows, plex S2 with last 10
  mutate(PlexID = paste0("S", rep(1:2, each = 10)),
         QuantBlock = 1) %>% 
  group_by(PlexID) %>% 
  # Within each of the two PlexID groups, create unique reporter aliases
  # and measurement names. ReporterAlias is "ref" for channel 131,
  # and MeasurementName is NA so it is not included in the cross-tab.
  mutate(ReporterAlias = paste(PlexID, 1:n(), sep = "_"),
         ReporterAlias = ifelse(ReporterName == "131", "ref", ReporterAlias),
         MeasurementName = ifelse(ReporterAlias == "ref", 
                                  NA, ReporterAlias)) %>% 
  ungroup() # stop grouping by PlexID
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:20em; "><table class="table table-hover table-condensed" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">(\#tab:samples-table)<left>Samples</left>
</caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> ReporterName </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> PlexID </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> QuantBlock </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> ReporterAlias </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> MeasurementName </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 126 </td>
   <td style="text-align:left;"> S1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> S1_1 </td>
   <td style="text-align:left;"> S1_1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 127N </td>
   <td style="text-align:left;"> S1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> S1_2 </td>
   <td style="text-align:left;"> S1_2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 127C </td>
   <td style="text-align:left;"> S1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> S1_3 </td>
   <td style="text-align:left;"> S1_3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 128N </td>
   <td style="text-align:left;"> S1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> S1_4 </td>
   <td style="text-align:left;"> S1_4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 128C </td>
   <td style="text-align:left;"> S1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> S1_5 </td>
   <td style="text-align:left;"> S1_5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 129N </td>
   <td style="text-align:left;"> S1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> S1_6 </td>
   <td style="text-align:left;"> S1_6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 129C </td>
   <td style="text-align:left;"> S1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> S1_7 </td>
   <td style="text-align:left;"> S1_7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 130N </td>
   <td style="text-align:left;"> S1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> S1_8 </td>
   <td style="text-align:left;"> S1_8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 130C </td>
   <td style="text-align:left;"> S1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> S1_9 </td>
   <td style="text-align:left;"> S1_9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 131 </td>
   <td style="text-align:left;"> S1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> ref </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 126 </td>
   <td style="text-align:left;"> S2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> S2_1 </td>
   <td style="text-align:left;"> S2_1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 127N </td>
   <td style="text-align:left;"> S2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> S2_2 </td>
   <td style="text-align:left;"> S2_2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 127C </td>
   <td style="text-align:left;"> S2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> S2_3 </td>
   <td style="text-align:left;"> S2_3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 128N </td>
   <td style="text-align:left;"> S2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> S2_4 </td>
   <td style="text-align:left;"> S2_4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 128C </td>
   <td style="text-align:left;"> S2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> S2_5 </td>
   <td style="text-align:left;"> S2_5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 129N </td>
   <td style="text-align:left;"> S2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> S2_6 </td>
   <td style="text-align:left;"> S2_6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 129C </td>
   <td style="text-align:left;"> S2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> S2_7 </td>
   <td style="text-align:left;"> S2_7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 130N </td>
   <td style="text-align:left;"> S2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> S2_8 </td>
   <td style="text-align:left;"> S2_8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 130C </td>
   <td style="text-align:left;"> S2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> S2_9 </td>
   <td style="text-align:left;"> S2_9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 131 </td>
   <td style="text-align:left;"> S2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> ref </td>
   <td style="text-align:left;"> NA </td>
  </tr>
</tbody>
</table></div>

</br>

Table \@ref(tab:samples-table) shows the `samples` table.


### References {#global-references}

The reference can be a certain channel, average of multiple channels, or 1 (no reference). The general form is an expression with `ReporterAlias` names as variables. It is evaluated for each `PlexID`/`QuantBlock` combination and applied to divide reporter ion intensities within corresponding `PlexID`/`QuantBlock`. A reference is used to convert raw intensities to relative intensities.


```r
# Create references table
references <- samples %>% 
  filter(ReporterAlias == "ref") %>% 
  # Select required columns and rename ReporterAlias to Reference
  select(PlexID, QuantBlock, Reference = ReporterAlias)
```

<table class="table table-hover table-condensed" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">(\#tab:references-table)<left>References</left>
</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> PlexID </th>
   <th style="text-align:right;"> QuantBlock </th>
   <th style="text-align:left;"> Reference </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> S1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> ref </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> ref </td>
  </tr>
</tbody>
</table>

</br>

Table \@ref(tab:references-table) shows the `references` table. The code to use the geometric average instead of a single channel as the reference is shown below. The geometric average is the product of the reporter ion channels to the power of (1/number of channels). For each `PlexID` group, collapse the vector of reporter ion names with `*`, surround them in parentheses, and raise to the power of (1/number of channels).

**Note:** if using the geometric average or no reference, make sure the samples table is correct. There should not be any `ReporterAlias` that are "ref" or `MeasurementName` that are `NA`.


```r
## Example of how to use the geometric average as reference - not run
references <- samples %>%
  group_by(PlexID, QuantBlock) %>%
  summarise(Reference = sprintf("(%s)^(1/%d)", 
                                paste(ReporterAlias, collapse = "*"), 
                                n()),
            .groups = "keep")
```


```r
## Example of how to set the reference to 1 - not run
references <- samples %>% 
  distinct(PlexID, QuantBlock) %>% 
  mutate(Reference = 1)
```

Now that we have the three study design tables, we should save them.


```r
# Save study design tables with write.table
write.table(fractions, file = "data/3442_fractions.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(samples, file = "data/3442_samples.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(references, file = "data/3442_references.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
```

Once the study design tables have been saved to text files, it is good practice to make them available to others. To do so, navigate to the Share Path provided in the DMS Data Package Detail Report (shown in Figure \@ref(fig:DMS-share-path)), and copy the three study design files to this location. This allows them to be accessed by others with the `get_study_design_by_dataset_package` function in the future.

<div class="figure" style="text-align: center">
<img src="images/data_package_share_path.PNG" alt="Location of the Share Path used to add the study design tables." width="75%" />
<p class="caption">(\#fig:DMS-share-path)Location of the Share Path used to add the study design tables.</p>
</div>


## Create Quantitative Cross-tab {#global-quant-crosstab}

This is the step where MS/MS IDs and reporter ions are linked together and aggregated to the peptide or accession (i.e. protein) level. To retain protein IDs while aggregating to peptide level, set `aggregation_level <- c("accession","peptide")`. The aggregation level can be any column or combination of columns in `psms(msnid)`. If specified by the study design tables, the intensities are converted to relative intensities by dividing by a reference. Then, they are log$_2$-transformed.


```r
# Create protein-level cross-tab by aggregating to accession level
crosstab <- create_crosstab(msnid = msnid, 
                            reporter_intensities = masic_data,
                            aggregation_level = "accession",
                            fractions = fractions, 
                            samples = samples, 
                            references = references)
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:90%; "><table class="table table-hover table-condensed" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">(\#tab:cross-tab)<left>First 6 rows of the cross-tab.</left>
</caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_3 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_4 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_5 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_6 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_7 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_8 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_9 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_1 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_2 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_3 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_4 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_5 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_6 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_7 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_8 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_9 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> AP_004893.1 </td>
   <td style="text-align:right;"> 0.1419768 </td>
   <td style="text-align:right;"> 0.1653552 </td>
   <td style="text-align:right;"> 0.7628195 </td>
   <td style="text-align:right;"> 0.9453172 </td>
   <td style="text-align:right;"> 0.8662554 </td>
   <td style="text-align:right;"> -1.9294467 </td>
   <td style="text-align:right;"> -0.6460065 </td>
   <td style="text-align:right;"> -1.2831873 </td>
   <td style="text-align:right;"> -0.4321433 </td>
   <td style="text-align:right;"> -1.0271227 </td>
   <td style="text-align:right;"> 0.4883309 </td>
   <td style="text-align:right;"> -0.9390945 </td>
   <td style="text-align:right;"> -0.7029685 </td>
   <td style="text-align:right;"> -1.7148628 </td>
   <td style="text-align:right;"> -0.1912097 </td>
   <td style="text-align:right;"> -0.8794712 </td>
   <td style="text-align:right;"> -0.2440478 </td>
   <td style="text-align:right;"> 0.3964607 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AP_004894.1 </td>
   <td style="text-align:right;"> 0.7947114 </td>
   <td style="text-align:right;"> -0.3151990 </td>
   <td style="text-align:right;"> -0.0913574 </td>
   <td style="text-align:right;"> 0.1974134 </td>
   <td style="text-align:right;"> 0.3033858 </td>
   <td style="text-align:right;"> -0.1750536 </td>
   <td style="text-align:right;"> -0.3527197 </td>
   <td style="text-align:right;"> -1.1762004 </td>
   <td style="text-align:right;"> -0.6438817 </td>
   <td style="text-align:right;"> -0.5124954 </td>
   <td style="text-align:right;"> -0.4428327 </td>
   <td style="text-align:right;"> -0.2364175 </td>
   <td style="text-align:right;"> -0.6711809 </td>
   <td style="text-align:right;"> -1.3730408 </td>
   <td style="text-align:right;"> -0.7462995 </td>
   <td style="text-align:right;"> -1.3515366 </td>
   <td style="text-align:right;"> -0.2227493 </td>
   <td style="text-align:right;"> -0.8338103 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AP_004895.1 </td>
   <td style="text-align:right;"> 0.2078433 </td>
   <td style="text-align:right;"> -0.6089756 </td>
   <td style="text-align:right;"> -0.2867209 </td>
   <td style="text-align:right;"> -0.3840271 </td>
   <td style="text-align:right;"> -0.1162062 </td>
   <td style="text-align:right;"> -0.6908468 </td>
   <td style="text-align:right;"> -1.1240967 </td>
   <td style="text-align:right;"> -0.7140383 </td>
   <td style="text-align:right;"> -0.6652575 </td>
   <td style="text-align:right;"> 0.2843676 </td>
   <td style="text-align:right;"> -0.1312555 </td>
   <td style="text-align:right;"> -0.1477038 </td>
   <td style="text-align:right;"> -0.4352950 </td>
   <td style="text-align:right;"> -0.6371609 </td>
   <td style="text-align:right;"> -0.6150788 </td>
   <td style="text-align:right;"> -0.6819180 </td>
   <td style="text-align:right;"> -0.1602120 </td>
   <td style="text-align:right;"> -0.3978979 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AP_004896.1 </td>
   <td style="text-align:right;"> -0.1494849 </td>
   <td style="text-align:right;"> -0.7314368 </td>
   <td style="text-align:right;"> -0.3664339 </td>
   <td style="text-align:right;"> -0.5352280 </td>
   <td style="text-align:right;"> -0.1742391 </td>
   <td style="text-align:right;"> -1.0372327 </td>
   <td style="text-align:right;"> -1.2945071 </td>
   <td style="text-align:right;"> -0.8299749 </td>
   <td style="text-align:right;"> -0.7060783 </td>
   <td style="text-align:right;"> 0.1939540 </td>
   <td style="text-align:right;"> -0.1688422 </td>
   <td style="text-align:right;"> -0.2274358 </td>
   <td style="text-align:right;"> -0.4222698 </td>
   <td style="text-align:right;"> -0.5251264 </td>
   <td style="text-align:right;"> -0.6741064 </td>
   <td style="text-align:right;"> -0.6543311 </td>
   <td style="text-align:right;"> -0.0441485 </td>
   <td style="text-align:right;"> -0.3994149 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AP_004898.1 </td>
   <td style="text-align:right;"> 0.0362964 </td>
   <td style="text-align:right;"> 0.4252227 </td>
   <td style="text-align:right;"> 0.7497227 </td>
   <td style="text-align:right;"> 1.1580326 </td>
   <td style="text-align:right;"> 0.4913660 </td>
   <td style="text-align:right;"> -0.3640632 </td>
   <td style="text-align:right;"> 0.1211536 </td>
   <td style="text-align:right;"> -0.8291744 </td>
   <td style="text-align:right;"> -0.3019505 </td>
   <td style="text-align:right;"> -0.8407749 </td>
   <td style="text-align:right;"> -0.4130732 </td>
   <td style="text-align:right;"> -0.2796091 </td>
   <td style="text-align:right;"> -0.9449498 </td>
   <td style="text-align:right;"> -1.5747761 </td>
   <td style="text-align:right;"> -0.1774225 </td>
   <td style="text-align:right;"> -1.8439756 </td>
   <td style="text-align:right;"> -0.4175363 </td>
   <td style="text-align:right;"> -1.1083199 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AP_004899.1 </td>
   <td style="text-align:right;"> 0.7140968 </td>
   <td style="text-align:right;"> -0.3732752 </td>
   <td style="text-align:right;"> -0.1781542 </td>
   <td style="text-align:right;"> -0.0615626 </td>
   <td style="text-align:right;"> 0.3494902 </td>
   <td style="text-align:right;"> -0.8550940 </td>
   <td style="text-align:right;"> -2.1679002 </td>
   <td style="text-align:right;"> -1.4519278 </td>
   <td style="text-align:right;"> -0.9026145 </td>
   <td style="text-align:right;"> -0.3158081 </td>
   <td style="text-align:right;"> -0.4644758 </td>
   <td style="text-align:right;"> -0.4056811 </td>
   <td style="text-align:right;"> -0.9023044 </td>
   <td style="text-align:right;"> -0.2805080 </td>
   <td style="text-align:right;"> -0.8052899 </td>
   <td style="text-align:right;"> -1.0482424 </td>
   <td style="text-align:right;"> -0.3959923 </td>
   <td style="text-align:right;"> -0.6675429 </td>
  </tr>
</tbody>
</table></div>

</br>

Now that we have the cross-tab, we should save it.


```r
# Save cross-tab
write.table(crosstab, file = "data/3442_global_crosstab.txt",
            sep = "\t", quote = FALSE, row.names = TRUE)
```

We will also save the proteins (row names) of this cross-tab in order to demonstrate prioritized inference later on.


```r
# Save global proteins
global_proteins <- rownames(crosstab)
save(global_proteins, file = "data/3442_global_proteins.RData")
```


## Create MSnSet {#global-msnset}

The `create_msnset` function can be used to easily create an MSnSet from the cross-tab and samples tables. More details about MSnSets will be added in a separate section at a later date. For now, read the documentation with `help("MSnSet")` or `?MSnSet`.


```r
# Create MSnSet
m <- create_msnset(crosstab = crosstab, samples = samples)
m
```

```
## MSnSet (storageMode: lockedEnvironment)
## assayData: 5173 features, 18 samples 
##   element names: exprs 
## protocolData: none
## phenoData
##   sampleNames: S1_1 S1_2 ... S2_9 (18 total)
##   varLabels: ReporterName PlexID ... MeasurementName (5 total)
##   varMetadata: labelDescription
## featureData: none
## experimentData: use 'experimentData(object)'
## Annotation:  
## - - - Processing information - - -
##  MSnbase version: 2.18.0
```



```r
# Save global MSnSet
save(m, file = "data/global_msnset.RData", compress = TRUE)
```



