# Isobaric Quantification Pipelines {#isobaric-quantification-pipelines}

<!---
TODO:
* Overview of isobaric labeling
--->

Section summary to be added later. Please continue to subsections.


## Global Proteomics Data {#global-proteomics-data}

This pipeline shows how to process TMT data that is processed outside of PNNL's DMS. Section \@ref(phosphoproteomics-data) shows how to process data from the DMS. For convenience, the results of MS-GF+ and MASIC processing are provided in a companion `PlexedPiperTestData` package. For this section, we need three packages: `PlexedPiper` for isobaric quantification, `PlexedPiperTestData`, and `dplyr` to manipulate data frames.

<!---
For these examples, we need two packages: PlexedPiper and PlexedPiperTestData. The former is used to interface with the DMS, while the latter is tailored to isobaric quantification.
--->





```r
# Setup
library(PlexedPiper)
library(PlexedPiperTestData)
library(dplyr)
```




### Prepare MS/MS Identifications {#prepare-MS2-IDs-global}

#### Read MS-GF+ Data {-}


```r
# Get file path
path_to_MSGF_results <- system.file("extdata/global/msgf_output", 
                                    package = "PlexedPiperTestData")
# Read MS-GF+ data from path
msnid <- read_msgf_data(path_to_MSGF_results)
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

`msnid` consists of 4 spectrum files (datasets), and contains a total of 1,156,754 peptide-spectrum-matches (PSMs), 511,617 total peptides, and 128,378 total accessions (proteins). The reported FDR is the empirical **false-discovery rate**, which is calculated as the ratio of the number of false (decoy) PSMs, peptides, or accessions to their true (non-decoy) counterparts. Calculation of these counts and their FDRs is shown below.


```r
# Calculating the counts and FDRs from the show() output ---
# Spectrum Files:
# Count
psms(msnid) %>% 
  distinct(Dataset) %>%
  nrow() # 48

# PSMs:
# Count
psms(msnid) %>% 
  distinct(Dataset, Scan, peptide, isDecoy) %>%
  # Assign intermediate to variable
  assign("x_psm", ., envir = globalenv()) %>% 
  nrow() # 1156754
# FDR
nrow(x_psm[x_psm$isDecoy == TRUE, ]) / 
  nrow(x_psm[x_psm$isDecoy == FALSE, ]) 
# 0.3127463 = 31%

# peptides:
# Count
psms(msnid) %>% 
  distinct(peptide, isDecoy) %>% 
  assign("x_peptide", ., envir = globalenv()) %>%
  nrow() # 511617
# FDR
nrow(x_peptide[x_peptide$isDecoy == TRUE, ]) / 
  nrow(x_peptide[x_peptide$isDecoy == FALSE, ])
# 0.611245 = 61%

# Accessions:
# Count
length(accessions(msnid)) # or
psms(msnid) %>% 
  distinct(accession, isDecoy) %>% 
  assign("x_acc", ., envir = globalenv()) %>%
  nrow() # 128378
# FDR
nrow(x_acc[x_acc$isDecoy == TRUE, ]) / 
  nrow(x_acc[x_acc$isDecoy == FALSE, ])
# 0.9827024 = 98%
```

Now that we have an `MSnID` object, we need to process it. 

#### Correct Isotope Selection Error {-}

Occasionally, the instrument selects a peak with +1 or more C13 atoms, rather than the monoisotopic (lowest mass) peak. While MS-FG+ is still capable of correctly identifying those, the downstream calculations of mass measurement error need to be fixed. The `correct_peak_selection` method corrects for these mass measurement errors.


```r
# Correct for isotope selection error
msnid <- correct_peak_selection(msnid)
```

#### Remove Contaminants {-}

Now, we will remove contaminants such as the pig trypsin that was used for protein digestion. We can use `grepl` to search for all accessions that contain "Contaminant".


```r
# All unique contaminants
unique(msnid$accession[grepl("Contaminant", msnid$accession)])
```

```
##  [1] "Contaminant_K2C1_HUMAN"      "Contaminant_K1C9_HUMAN"     
##  [3] "Contaminant_ALBU_HUMAN"      "Contaminant_ALBU_BOVIN"     
##  [5] "Contaminant_TRYP_PIG"        "Contaminant_K1C10_HUMAN"    
##  [7] "XXX_Contaminant_K1C9_HUMAN"  "Contaminant_K22E_HUMAN"     
##  [9] "Contaminant_Trypa3"          "Contaminant_Trypa5"         
## [11] "XXX_Contaminant_K1C10_HUMAN" "XXX_Contaminant_K22E_HUMAN" 
## [13] "XXX_Contaminant_K2C1_HUMAN"  "Contaminant_TRYP_BOVIN"     
## [15] "XXX_Contaminant_ALBU_HUMAN"  "XXX_Contaminant_ALBU_BOVIN" 
## [17] "XXX_Contaminant_TRYP_BOVIN"  "Contaminant_CTRB_BOVIN"     
## [19] "Contaminant_Trypa1"          "Contaminant_Trypa6"         
## [21] "Contaminant_CTRA_BOVIN"      "XXX_Contaminant_TRYP_PIG"   
## [23] "XXX_Contaminant_CTRB_BOVIN"  "XXX_Contaminant_CTRA_BOVIN" 
## [25] "Contaminant_Trypa2"
```

To remove contaminants, we use `apply_filter` with an appropriate character string that tells the function what rows to keep. In this case, we keep rows where the accession does not contain "Contaminant". We will use `show` to see how the counts change.


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

#### MS/MS ID Filter: Peptide Level {-}

The next step is to use the `PepQValue` column from the MSnID object and the absolute deviation of the mass measurement error of parent ions (in ppm) to maximize the number of PSMs while ensuring that the empirical peptide-level FDR is at most 1%.


```r
# 1% FDR filter at the peptide level
msnid <- filter_msgf_data(msnid,
                          level = "peptide",
                          fdr.max = 0.01)
show(msnid)
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  48 
## #PSMs: 464542 at 0.45 % FDR
## #peptides: 96493 at 1 % FDR
## #accessions: 27120 at 9.2 % FDR
```

We can see that filtering drastically reduces the number of PSMs, and the empirical peptide-level FDR is now 1%. However, notice that the empirical protein-level FDR is still fairly high.

#### MS/MS ID Filter: Protein Level {-}

A while ago, the proteomics field established the hard-and-fast two-peptides-per-protein rule. That is, we can not be confident if a protein is identified by the detection of only one peptide. This rule penalizes short proteins and doesn't consider that there are some very long proteins (e.g. Titin 3.8 MDa) that easily have more then two matching peptides even in reversed sequence. Thus, we propose to normalize the number of peptides per protein length and use that as a filtering criterion.

We need the FASTA (pronounced FAST-AYE) file to get the length of each protein, which we can then use to calculate the associated number of peptides per 1000 amino acids. This new `peptides_per_1000aa` column is used to filter the MSnID object so that the empirical accession-level FDR is at most 1%.


```r
# Get path to FASTA file
path_to_FASTA <- system.file(
  "extdata/Rattus_norvegicus_NCBI_RefSeq_2018-04-10.fasta.gz", 
  package = "PlexedPiperTestData"
)

# Compute number of peptides per 1000 amino acids
msnid <- compute_num_peptides_per_1000aa(msnid, path_to_FASTA)

# 1% FDR filter at the protein level
msnid <- filter_msgf_data(msnid,
                          level = "accession",
                          fdr.max = 0.01)
show(msnid)
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  48 
## #PSMs: 458097 at 0.16 % FDR
## #peptides: 92024 at 0.32 % FDR
## #accessions: 15620 at 0.98 % FDR
```

#### Inference of Parsimonious Protein Set {-}

The situation when a certain peptide sequence matches multiple proteins adds complication to the downstream quantitative analysis, as it is not clear which protein this peptide is originating from. There are common ways for dealing with this. One is to simply retain uniquely matching peptides and discard shared peptides (`unique_only = TRUE`). Alternatively (in case of `unique_only = FALSE`) assign the shared peptides to the proteins with the larger number of uniquely mapping peptides. If there is a choice between multiple proteins with equal numbers of uniquely mapping peptides, the shared peptides are assigned to the first protein according to alphanumeric order. This step could be done prior to filtering at the accession level, but the removal of an accession will completely remove its associated peptides.

<!---
TODO:
* Include a visual for inference of the parsimonious protein set.
--->


```r
# Inference of parsimonious protein set
msnid <- infer_parsimonious_accessions(msnid, unique_only = FALSE)
show(msnid)
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  48 
## #PSMs: 445003 at 0.15 % FDR
## #peptides: 90466 at 0.27 % FDR
## #accessions: 5246 at 1.1 % FDR
```

Notice that the protein-level FDR increased above the acceptable threshold, so we need to reapply the filter.


```r
# 1% FDR filter at the protein level
msnid <- filter_msgf_data(msnid,
                          level = "accession",
                          fdr.max = 0.01)
show(msnid)
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  48 
## #PSMs: 444857 at 0.14 % FDR
## #peptides: 90357 at 0.26 % FDR
## #accessions: 5211 at 0.99 % FDR
```

Once all filtering is done, we can remove the decoy accessions. We use the `apply_filter` function again and only keep entries where `isDecoy` is `FALSE`.


```r
# Remove Decoy Accessions
msnid <- apply_filter(msnid, "!isDecoy")
show(msnid)
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  48 
## #PSMs: 444216 at 0 % FDR
## #peptides: 90126 at 0 % FDR
## #accessions: 5160 at 0 % FDR
```

After processing, we are left with 318,448 PSMs, 81,048 peptides, and 5,143 proteins. The empirical FDRs are the same as before, but can not be calculated because we removed the decoys.

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:20em; overflow-x: scroll; width:100%; "><table class="table" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">(\#tab:global-msnid-table)<left>First 10 rows of the processed MS-GF+ results.</left>
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
   <td style="text-align:right;"> 25.769 </td>
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
   <td style="text-align:right;"> 25.769 </td>
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
   <td style="text-align:right;"> 25.769 </td>
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
   <td style="text-align:right;"> 25.769 </td>
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
   <td style="text-align:right;"> 580.844 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_08_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 18244 </td>
   <td style="text-align:right;"> 10302 </td>
   <td style="text-align:left;"> HCD </td>
   <td style="text-align:right;"> 3724 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 653.695 </td>
   <td style="text-align:right;"> -0.003 </td>
   <td style="text-align:right;"> -1.649 </td>
   <td style="text-align:right;"> 1958.071 </td>
   <td style="text-align:left;"> A.AAAAATEQQGSNGPVK.K </td>
   <td style="text-align:left;"> NP_001177997.1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 106 </td>
   <td style="text-align:right;"> 85 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NP_001177997.1 </td>
   <td style="text-align:right;"> 653.362 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 653.361 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_08_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 10302 </td>
   <td style="text-align:left;"> AAAAATEQQGSNGPVK </td>
   <td style="text-align:right;"> 3.432 </td>
   <td style="text-align:right;"> 1.685 </td>
   <td style="text-align:right;"> 62.500 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_16_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 13019 </td>
   <td style="text-align:right;"> 8116 </td>
   <td style="text-align:left;"> HCD </td>
   <td style="text-align:right;"> 3103 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 613.301 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.100 </td>
   <td style="text-align:right;"> 1837.888 </td>
   <td style="text-align:left;"> R.AAAADGEPLHNEEER.T </td>
   <td style="text-align:left;"> NP_001099982.1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 112 </td>
   <td style="text-align:right;"> 108 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> NP_001099982.1 </td>
   <td style="text-align:right;"> 613.301 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 613.301 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_16_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 8116 </td>
   <td style="text-align:left;"> AAAADGEPLHNEEER </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> 0.105 </td>
   <td style="text-align:right;"> 14.749 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_16_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 18313 </td>
   <td style="text-align:right;"> 8217 </td>
   <td style="text-align:left;"> HCD </td>
   <td style="text-align:right;"> 3907 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 613.302 </td>
   <td style="text-align:right;"> 0.003 </td>
   <td style="text-align:right;"> 1.692 </td>
   <td style="text-align:right;"> 1837.888 </td>
   <td style="text-align:left;"> R.AAAADGEPLHNEEER.T </td>
   <td style="text-align:left;"> NP_001099982.1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 90 </td>
   <td style="text-align:right;"> 68 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.001 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> NP_001099982.1 </td>
   <td style="text-align:right;"> 613.301 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 613.302 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S2_16_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 8217 </td>
   <td style="text-align:left;"> AAAADGEPLHNEEER </td>
   <td style="text-align:right;"> 3.000 </td>
   <td style="text-align:right;"> 1.703 </td>
   <td style="text-align:right;"> 14.749 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_24_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 16441 </td>
   <td style="text-align:right;"> 6833 </td>
   <td style="text-align:left;"> HCD </td>
   <td style="text-align:right;"> 3438 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 795.928 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 1590.849 </td>
   <td style="text-align:left;"> A.AAAAEAESGGGGGK.K </td>
   <td style="text-align:left;"> NP_001128630.1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 230 </td>
   <td style="text-align:right;"> 192 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> NP_001128630.1 </td>
   <td style="text-align:right;"> 795.928 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 795.928 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_24_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 6833 </td>
   <td style="text-align:left;"> AAAAEAESGGGGGK </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> 0.104 </td>
   <td style="text-align:right;"> 580.844 </td>
  </tr>
</tbody>
</table></div>

</br>

Table \@ref(tab:global-msnid-table) shows the first 10 rows of the processed MS-GF+ output.


### Prepare Reporter Ion Intensities {#reporter-ion-intensities}

#### Read MASIC Output {-}

MASIC is a tool for extracting ion intensities. With proper parameter settings, it can be used for extracting TMT (or iTRAQ) reporter ion intensities. In addition, it reports a number of other helpful metrics. Notably, the interference score at the parent ion level and the signal-to-noise ratio (S/N) at the reporter ion level (computed by Thermo software). The interference score reflects the proportion of the ion population that was isolated for fragmentation that is due to the targeted ion. In other words, 1 - InterferenceScore is due to co-isolated species that have similar elution time and parent ion m/z.


```r
# Path to MASIC data
path_to_MASIC_results <- system.file("extdata/global/masic_output", 
                                     package = "PlexedPiperTestData")

# Read MASIC data
masic_data <- read_masic_data(path_to_MASIC_results, interference_score = TRUE)
```

Normally, this would display two progress bars in the console as the data is being fetched. However, the output was suppressed to save space.

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:20em; overflow-x: scroll; width:100%; "><table class="table" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">(\#tab:global-masic-table)<left>First 10 rows of the MASIC data.</left>
</caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Dataset </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ScanNumber </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Collision.Mode </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ParentIonMZ </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> BasePeakIntensity </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> BasePeakMZ </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ReporterIonIntensityMax </th>
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
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Weighted.Avg.Pct.Intensity.Correction </th>
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
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_126.128_Resolution </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_127.125_Resolution </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_127.131_Resolution </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_128.128_Resolution </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_128.134_Resolution </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_129.131_Resolution </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_129.138_Resolution </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_130.135_Resolution </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_130.141_Resolution </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Ion_131.138_Resolution </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ParentIonIndex </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> MZ </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> SurveyScanNumber </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> OptimalPeakApexScanNumber </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PeakApexOverrideParentIonIndex </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> CustomSICPeak </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PeakScanStart </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PeakScanEnd </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PeakScanMaxIntensity </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PeakMaxIntensity </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PeakSignalToNoiseRatio </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> FWHMInScans </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PeakArea </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ParentIonIntensity </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PeakBaselineNoiseLevel </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PeakBaselineNoiseStDev </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PeakBaselinePointsUsed </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> StatMomentsArea </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> CenterOfMassScan </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PeakStDev </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PeakSkew </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> PeakKSStat </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> StatMomentsDataCountUsed </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> InterferenceScore </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> hcd </td>
   <td style="text-align:right;"> 407.74 </td>
   <td style="text-align:right;"> 227695.44 </td>
   <td style="text-align:right;"> 407.741 </td>
   <td style="text-align:right;"> 92236.87 </td>
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
   <td style="text-align:right;"> 0 </td>
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
   <td style="text-align:right;"> 44102 </td>
   <td style="text-align:right;"> 42700 </td>
   <td style="text-align:right;"> 42100 </td>
   <td style="text-align:right;"> 41800 </td>
   <td style="text-align:right;"> 44404 </td>
   <td style="text-align:right;"> 40500 </td>
   <td style="text-align:right;"> 39500 </td>
   <td style="text-align:right;"> 36800 </td>
   <td style="text-align:right;"> 41100 </td>
   <td style="text-align:right;"> 42302 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 407.742 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> -1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 2901600 </td>
   <td style="text-align:right;"> 211.000 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 50422000 </td>
   <td style="text-align:right;"> 2579600 </td>
   <td style="text-align:right;"> 13750 </td>
   <td style="text-align:right;"> 97562 </td>
   <td style="text-align:right;"> 10113 </td>
   <td style="text-align:right;"> 47031000 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 5.68 </td>
   <td style="text-align:right;"> -0.189 </td>
   <td style="text-align:right;"> 0.353 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0.996 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> hcd </td>
   <td style="text-align:right;"> 396.72 </td>
   <td style="text-align:right;"> 59127.97 </td>
   <td style="text-align:right;"> 529.294 </td>
   <td style="text-align:right;"> 34294.90 </td>
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
   <td style="text-align:right;"> 0 </td>
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
   <td style="text-align:right;"> 42702 </td>
   <td style="text-align:right;"> 41100 </td>
   <td style="text-align:right;"> 37000 </td>
   <td style="text-align:right;"> 40400 </td>
   <td style="text-align:right;"> 43404 </td>
   <td style="text-align:right;"> 36400 </td>
   <td style="text-align:right;"> 39700 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 29800 </td>
   <td style="text-align:right;"> 41802 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 396.718 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> -1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 2181900 </td>
   <td style="text-align:right;"> 19.690 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 34508000 </td>
   <td style="text-align:right;"> 1690600 </td>
   <td style="text-align:right;"> 110841 </td>
   <td style="text-align:right;"> 1120000 </td>
   <td style="text-align:right;"> 10166 </td>
   <td style="text-align:right;"> 31578000 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 5.59 </td>
   <td style="text-align:right;"> -0.217 </td>
   <td style="text-align:right;"> 0.347 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:right;"> 0.993 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> hcd </td>
   <td style="text-align:right;"> 438.23 </td>
   <td style="text-align:right;"> 110444.82 </td>
   <td style="text-align:right;"> 362.224 </td>
   <td style="text-align:right;"> 14053.40 </td>
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
   <td style="text-align:right;"> 0 </td>
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
   <td style="text-align:right;"> 42006 </td>
   <td style="text-align:right;"> 40702 </td>
   <td style="text-align:right;"> 41402 </td>
   <td style="text-align:right;"> 40700 </td>
   <td style="text-align:right;"> 40400 </td>
   <td style="text-align:right;"> 38800 </td>
   <td style="text-align:right;"> 40200 </td>
   <td style="text-align:right;"> 38900 </td>
   <td style="text-align:right;"> 40400 </td>
   <td style="text-align:right;"> 41002 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 438.227 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 131 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 114 </td>
   <td style="text-align:right;"> 107 </td>
   <td style="text-align:right;"> 8255600 </td>
   <td style="text-align:right;"> 9.465 </td>
   <td style="text-align:right;"> 44 </td>
   <td style="text-align:right;"> 444610000 </td>
   <td style="text-align:right;"> 658727 </td>
   <td style="text-align:right;"> 872195 </td>
   <td style="text-align:right;"> 2620000 </td>
   <td style="text-align:right;"> 10129 </td>
   <td style="text-align:right;"> 343470000 </td>
   <td style="text-align:right;"> 82 </td>
   <td style="text-align:right;"> 22.80 </td>
   <td style="text-align:right;"> -0.626 </td>
   <td style="text-align:right;"> 1.099 </td>
   <td style="text-align:right;"> 89 </td>
   <td style="text-align:right;"> 1.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:left;"> hcd </td>
   <td style="text-align:right;"> 481.50 </td>
   <td style="text-align:right;"> 37082.72 </td>
   <td style="text-align:right;"> 206.466 </td>
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
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0 </td>
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
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 481.505 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> -1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 52 </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 401824 </td>
   <td style="text-align:right;"> 27.990 </td>
   <td style="text-align:right;"> 23 </td>
   <td style="text-align:right;"> 16244000 </td>
   <td style="text-align:right;"> 344491 </td>
   <td style="text-align:right;"> 14356 </td>
   <td style="text-align:right;"> 65777 </td>
   <td style="text-align:right;"> 10109 </td>
   <td style="text-align:right;"> 14899000 </td>
   <td style="text-align:right;"> 26 </td>
   <td style="text-align:right;"> 14.10 </td>
   <td style="text-align:right;"> 0.054 </td>
   <td style="text-align:right;"> 0.504 </td>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 1.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:left;"> hcd </td>
   <td style="text-align:right;"> 549.28 </td>
   <td style="text-align:right;"> 21077.05 </td>
   <td style="text-align:right;"> 128.129 </td>
   <td style="text-align:right;"> 21077.05 </td>
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
   <td style="text-align:right;"> 0 </td>
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
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 40302 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 42102 </td>
   <td style="text-align:right;"> 46600 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 40300 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 549.279 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> -1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 363656 </td>
   <td style="text-align:right;"> 0.700 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 5941200 </td>
   <td style="text-align:right;"> 347071 </td>
   <td style="text-align:right;"> 519640 </td>
   <td style="text-align:right;"> 1990000 </td>
   <td style="text-align:right;"> 10109 </td>
   <td style="text-align:right;"> 2528900 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 2.49 </td>
   <td style="text-align:right;"> -0.583 </td>
   <td style="text-align:right;"> 0.430 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 1.000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:left;"> hcd </td>
   <td style="text-align:right;"> 388.72 </td>
   <td style="text-align:right;"> 40605.85 </td>
   <td style="text-align:right;"> 356.719 </td>
   <td style="text-align:right;"> 8087.76 </td>
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
   <td style="text-align:right;"> 0 </td>
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
   <td style="text-align:right;"> 40000 </td>
   <td style="text-align:right;"> 26400 </td>
   <td style="text-align:right;"> 30400 </td>
   <td style="text-align:right;"> 40400 </td>
   <td style="text-align:right;"> 44300 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 28800 </td>
   <td style="text-align:right;"> 28500 </td>
   <td style="text-align:right;"> 38700 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 388.720 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> -1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 478135 </td>
   <td style="text-align:right;"> 28.710 </td>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:right;"> 10718000 </td>
   <td style="text-align:right;"> 291189 </td>
   <td style="text-align:right;"> 16653 </td>
   <td style="text-align:right;"> 142562 </td>
   <td style="text-align:right;"> 10142 </td>
   <td style="text-align:right;"> 9961100 </td>
   <td style="text-align:right;"> 16 </td>
   <td style="text-align:right;"> 8.60 </td>
   <td style="text-align:right;"> -0.051 </td>
   <td style="text-align:right;"> 0.283 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 0.969 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:left;"> hcd </td>
   <td style="text-align:right;"> 403.21 </td>
   <td style="text-align:right;"> 40667.64 </td>
   <td style="text-align:right;"> 403.251 </td>
   <td style="text-align:right;"> 5860.96 </td>
   <td style="text-align:right;"> 4991.79 </td>
   <td style="text-align:right;"> 1274.12 </td>
   <td style="text-align:right;"> 5860.96 </td>
   <td style="text-align:right;"> 4699.99 </td>
   <td style="text-align:right;"> 4906.33 </td>
   <td style="text-align:right;"> 1782.18 </td>
   <td style="text-align:right;"> 4580.83 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 5.29 </td>
   <td style="text-align:right;"> 1.35 </td>
   <td style="text-align:right;"> 6.21 </td>
   <td style="text-align:right;"> 4.98 </td>
   <td style="text-align:right;"> 5.19 </td>
   <td style="text-align:right;"> 1.89 </td>
   <td style="text-align:right;"> 4.85 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 41600 </td>
   <td style="text-align:right;"> 21900 </td>
   <td style="text-align:right;"> 38604 </td>
   <td style="text-align:right;"> 39000 </td>
   <td style="text-align:right;"> 36800 </td>
   <td style="text-align:right;"> 26700 </td>
   <td style="text-align:right;"> 39200 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 403.214 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> -1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 54 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 309628 </td>
   <td style="text-align:right;"> 18.930 </td>
   <td style="text-align:right;"> 25 </td>
   <td style="text-align:right;"> 13671000 </td>
   <td style="text-align:right;"> 232858 </td>
   <td style="text-align:right;"> 16353 </td>
   <td style="text-align:right;"> 50315 </td>
   <td style="text-align:right;"> 3379 </td>
   <td style="text-align:right;"> 12441000 </td>
   <td style="text-align:right;"> 27 </td>
   <td style="text-align:right;"> 14.90 </td>
   <td style="text-align:right;"> 0.086 </td>
   <td style="text-align:right;"> 0.489 </td>
   <td style="text-align:right;"> 41 </td>
   <td style="text-align:right;"> 0.618 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 10 </td>
   <td style="text-align:left;"> hcd </td>
   <td style="text-align:right;"> 476.26 </td>
   <td style="text-align:right;"> 14753.41 </td>
   <td style="text-align:right;"> 577.352 </td>
   <td style="text-align:right;"> 6170.43 </td>
   <td style="text-align:right;"> 2622.09 </td>
   <td style="text-align:right;"> 1825.23 </td>
   <td style="text-align:right;"> 2987.08 </td>
   <td style="text-align:right;"> 1852.17 </td>
   <td style="text-align:right;"> 2309.64 </td>
   <td style="text-align:right;"> 1348.42 </td>
   <td style="text-align:right;"> 6170.43 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 1735.76 </td>
   <td style="text-align:right;"> 5115.81 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 2.61 </td>
   <td style="text-align:right;"> 1.80 </td>
   <td style="text-align:right;"> 2.95 </td>
   <td style="text-align:right;"> 1.81 </td>
   <td style="text-align:right;"> 2.26 </td>
   <td style="text-align:right;"> 1.31 </td>
   <td style="text-align:right;"> 5.99 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 1.67 </td>
   <td style="text-align:right;"> 4.93 </td>
   <td style="text-align:right;"> 30500 </td>
   <td style="text-align:right;"> 28200 </td>
   <td style="text-align:right;"> 29400 </td>
   <td style="text-align:right;"> 28800 </td>
   <td style="text-align:right;"> 29800 </td>
   <td style="text-align:right;"> 23300 </td>
   <td style="text-align:right;"> 39004 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 31200 </td>
   <td style="text-align:right;"> 39200 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 476.262 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> -1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 18 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 276349 </td>
   <td style="text-align:right;"> 3.873 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 4751000 </td>
   <td style="text-align:right;"> 277210 </td>
   <td style="text-align:right;"> 71355 </td>
   <td style="text-align:right;"> 309842 </td>
   <td style="text-align:right;"> 3361 </td>
   <td style="text-align:right;"> 3027500 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 5.27 </td>
   <td style="text-align:right;"> -0.363 </td>
   <td style="text-align:right;"> 0.376 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 0.616 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:left;"> hcd </td>
   <td style="text-align:right;"> 455.24 </td>
   <td style="text-align:right;"> 22331.15 </td>
   <td style="text-align:right;"> 207.407 </td>
   <td style="text-align:right;"> 1521.09 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 1521.09 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 1.54 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 31600 </td>
   <td style="text-align:right;"> 8 </td>
   <td style="text-align:right;"> 455.243 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> -1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 330250 </td>
   <td style="text-align:right;"> 0.043 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 4566900 </td>
   <td style="text-align:right;"> 330250 </td>
   <td style="text-align:right;"> 7747800 </td>
   <td style="text-align:right;"> 56200000 </td>
   <td style="text-align:right;"> 10118 </td>
   <td style="text-align:right;"> 2203700 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 2.99 </td>
   <td style="text-align:right;"> -0.545 </td>
   <td style="text-align:right;"> 0.373 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0.948 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:left;"> hcd </td>
   <td style="text-align:right;"> 411.72 </td>
   <td style="text-align:right;"> 15172.97 </td>
   <td style="text-align:right;"> 126.128 </td>
   <td style="text-align:right;"> 15172.97 </td>
   <td style="text-align:right;"> 15172.97 </td>
   <td style="text-align:right;"> 7654.57 </td>
   <td style="text-align:right;"> 10487.60 </td>
   <td style="text-align:right;"> 9163.71 </td>
   <td style="text-align:right;"> 9064.43 </td>
   <td style="text-align:right;"> 5550.98 </td>
   <td style="text-align:right;"> 6100.78 </td>
   <td style="text-align:right;"> 6920.02 </td>
   <td style="text-align:right;"> 5498.97 </td>
   <td style="text-align:right;"> 9974.23 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 16.25 </td>
   <td style="text-align:right;"> 8.21 </td>
   <td style="text-align:right;"> 11.25 </td>
   <td style="text-align:right;"> 9.85 </td>
   <td style="text-align:right;"> 9.75 </td>
   <td style="text-align:right;"> 5.98 </td>
   <td style="text-align:right;"> 6.57 </td>
   <td style="text-align:right;"> 7.47 </td>
   <td style="text-align:right;"> 5.93 </td>
   <td style="text-align:right;"> 10.83 </td>
   <td style="text-align:right;"> 41906 </td>
   <td style="text-align:right;"> 38300 </td>
   <td style="text-align:right;"> 41200 </td>
   <td style="text-align:right;"> 39900 </td>
   <td style="text-align:right;"> 40500 </td>
   <td style="text-align:right;"> 38800 </td>
   <td style="text-align:right;"> 39600 </td>
   <td style="text-align:right;"> 38600 </td>
   <td style="text-align:right;"> 41800 </td>
   <td style="text-align:right;"> 39002 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 411.719 </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> -1 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 344633 </td>
   <td style="text-align:right;"> 0.004 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;"> 10577000 </td>
   <td style="text-align:right;"> 257141 </td>
   <td style="text-align:right;"> 86211000 </td>
   <td style="text-align:right;"> 507000000 </td>
   <td style="text-align:right;"> 9997 </td>
   <td style="text-align:right;"> 1705500 </td>
   <td style="text-align:right;"> 34 </td>
   <td style="text-align:right;"> 1.93 </td>
   <td style="text-align:right;"> -0.722 </td>
   <td style="text-align:right;"> 0.397 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 1.000 </td>
  </tr>
</tbody>
</table></div>

</br>

Table \@ref(tab:global-masic-table) shows the first 10 rows of `masic_data`.


#### Filter MASIC Data {-}

Currently, we recommend keeping entries where at least 50% of the ion population is due to the targeted ion (interference score $\geq$ 0.5) and not filtering by S/N.


```r
# Filter MASIC data
masic_data <- filter_masic_data(masic_data, 
                                interference_score_threshold = 0.5,
                                s2n_threshold = 0)
```

### Create Study Design Tables {#fetch-study-design-tables}

To convert from PSMs and reporter ion intensities to meaningful quantitative data, it is necessary to know what are the samples in the reporter channels and what is the intended reference channel (or combination of channels). The entire study design is captured by three tables - fractions, samples, references. With newly processed data, these typically do not exist, and must be created. If the tables already exist, the code to access them is as follows.


```r
# Read tables from folder:
fractions <- read_tsv(system.file("extdata/study_design/fractions.txt", 
                                  package = "PlexedPiperTestData"))
samples <- read_tsv(system.file("extdata/study_design/samples.txt", 
                                package = "PlexedPiperTestData"))
references <- read_tsv(system.file("extdata/study_design/references.txt", 
                                   package = "PlexedPiperTestData"))

# If using a data package from the DMS:
study_design <- read_study_design_from_DMS(data_package_num)
fractions <- study_design$fractions
samples <- study_design$samples
references <- study_design$references
```


#### Fractions {-}

The fractions table consists of two columns: `Dataset` and `PlexID`. The `Dataset` column contains all of the unique datasets from `msnid$Dataset` or `masic_data$Dataset`. The `PlexID` column contains the plex ID associated with each dataset, and is typically an "S" followed by a number ("S1", "S2", etc.). We can extract the plex ID from the datasets. In this case, the plex ID always comes after "_W_", so we can use a regular expression (regex) to capture it (the first argument of `gsub`). The regex below says to capture an "S" followed by a single digit that appears after "_W_" and before an underscore.


```r
# Create fractions table
fractions <- data.frame(Dataset = unique(masic_data$Dataset)) %>% 
  mutate(PlexID = gsub(".*_W_(S\\d{1})_.*", "\\1", Dataset))
```

<table class="table" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">(\#tab:fractions-table)<left>First 10 rows of the fractions table.</left>
</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Dataset </th>
   <th style="text-align:left;"> PlexID </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_01_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_02_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_03_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_04_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_05_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_06_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_07_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_08_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_09_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_W_S1_10_12Oct17_Elm_AQ-17-09-02 </td>
   <td style="text-align:left;"> S1 </td>
  </tr>
</tbody>
</table>

</br>

Table \@ref(tab:fractions-table) shows the first 10 rows of `fractions`.


#### Samples {-}

The samples table contains columns `PlexID`, `QuantBlock`, `ReporterName`, `ReporterAlias`, and `MeasurementName`. The plex ID must be the same as the plex ID in the `fractions` table. `ReporterName` is the reporter ion name ("126", "127N", "127C", etc.). `ReporterAlias` is the intermediate between `ReporterName` and `MeasurementName` and is used for defining the reference. `MeasurementName` determines the column names for the final cross-tab, and must be unique and begin with a letter. `MeasurementName` is easily constructed by prepending `PlexID` to the `ReporterName`. Finally, `QuantBlock` can be thought of as a way of defining sub-plex. In a typical TMT experiment, `QuantBlock` is always 1. In case of 5 pairwise comparisons within TMT10, there will be 5 QuantBlocks (1-5) with a reference for each `QuantBlock`.

For this experiment, channel 131 will serve as the reference, so we set `MeasurementName` to `NA` when `ReporterName` is `"131"`. This will make the reference channel absent from the quantitative cross-tab. In cases where reporter ion intensities are not normalized by a reference channel (reference = 1) or they are normalized by the average of select channels, do not set any `MeasurementName` to `NA`.


```r
# TMT10 Reporter Converter table from MSnID package
conv <- reporter_converter$tmt10
plexes <- unique(fractions$PlexID)

# Reference channel
ref_channel <- "131"

# Create samples table
samples <- data.frame(PlexID = rep(plexes, each = nrow(conv)),  
                      ReporterName = rep(conv$ReporterName, 
                                         length(plexes))) %>% 
  mutate(ReporterAlias = sprintf("%s_%s", PlexID, ReporterName),
         MeasurementName = ReporterAlias,
         QuantBlock = 1,
         # Comment out this next part if the reference
         # is not one of the reporter ion channels.
         MeasurementName = ifelse(ReporterName == ref_channel,
                                  NA, MeasurementName)
  )
```

<table class="table" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">(\#tab:samples-table)<left>First 10 rows of the samples table.</left>
</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> PlexID </th>
   <th style="text-align:left;"> ReporterName </th>
   <th style="text-align:left;"> ReporterAlias </th>
   <th style="text-align:left;"> MeasurementName </th>
   <th style="text-align:right;"> QuantBlock </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> S1 </td>
   <td style="text-align:left;"> 126 </td>
   <td style="text-align:left;"> S1_126 </td>
   <td style="text-align:left;"> S1_126 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S1 </td>
   <td style="text-align:left;"> 127N </td>
   <td style="text-align:left;"> S1_127N </td>
   <td style="text-align:left;"> S1_127N </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S1 </td>
   <td style="text-align:left;"> 127C </td>
   <td style="text-align:left;"> S1_127C </td>
   <td style="text-align:left;"> S1_127C </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S1 </td>
   <td style="text-align:left;"> 128N </td>
   <td style="text-align:left;"> S1_128N </td>
   <td style="text-align:left;"> S1_128N </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S1 </td>
   <td style="text-align:left;"> 128C </td>
   <td style="text-align:left;"> S1_128C </td>
   <td style="text-align:left;"> S1_128C </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S1 </td>
   <td style="text-align:left;"> 129N </td>
   <td style="text-align:left;"> S1_129N </td>
   <td style="text-align:left;"> S1_129N </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S1 </td>
   <td style="text-align:left;"> 129C </td>
   <td style="text-align:left;"> S1_129C </td>
   <td style="text-align:left;"> S1_129C </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S1 </td>
   <td style="text-align:left;"> 130N </td>
   <td style="text-align:left;"> S1_130N </td>
   <td style="text-align:left;"> S1_130N </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S1 </td>
   <td style="text-align:left;"> 130C </td>
   <td style="text-align:left;"> S1_130C </td>
   <td style="text-align:left;"> S1_130C </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S1 </td>
   <td style="text-align:left;"> 131 </td>
   <td style="text-align:left;"> S1_131 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
</tbody>
</table>

</br>

Table \@ref(tab:samples-table) shows the first 10 rows of `samples`.


#### References {-}

Reference can be a certain channel, average of multiple channels, or 1. The general form is an expression with `ReporterAlias` names as variables. It is evaluated for each `PlexID`/`QuantBlock` combination and applied to divide reporter ion intensities within corresponding `PlexID`/`QuantBlock`.


```r
# Create references table
references <- samples %>% 
  # Filter to reference channel
  filter(ReporterName == ref_channel) %>% 
  # Select required columns and rename ReporterAlias to Reference
  select(PlexID, Reference = ReporterAlias, QuantBlock)
```

<table class="table" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">(\#tab:references-table)<left>References table.</left>
</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> PlexID </th>
   <th style="text-align:left;"> Reference </th>
   <th style="text-align:right;"> QuantBlock </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> S1 </td>
   <td style="text-align:left;"> S1_131 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> S2 </td>
   <td style="text-align:left;"> S2_131 </td>
   <td style="text-align:right;"> 1 </td>
  </tr>
</tbody>
</table>

</br>

Table \@ref(tab:references-table) shows the first 10 rows of `references`. The code to use the geometric average instead of a single channel as the reference is shown below. The geometric average is the product of the reporter ion channels to the power of (1/number of channels). For each `PlexID` group, collapse the vector of reporter ion names with `*`, surround them in parentheses, and raise to the power of (1/number of channels).


```r
# Use geometric average as reference
references <- samples %>%
  group_by(PlexID, QuantBlock) %>%
  summarise(Reference = sprintf("(%s)^(1/%d)", 
                                paste(ReporterAlias, collapse = "*"), n()))

# Do not normalize by reference channel (use 1 as the reference)
references <- samples %>% 
  distinct(PlexID, QuantBlock) %>% 
  mutate(Reference = 1)
```

Now that we have the three study design tables, we should save them.

<!---
TODO:
Explain how to add the study design tables to the DMS.
--->


```r
# Save study design tables
write.table(fractions, file = "fractions.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(samples, file = "samples.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(references, file = "references.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
```


### Create Quantitative Cross-tab {#global-quant-crosstab}

This is the final step where MS/MS IDs and reporter ions are linked together and aggregated to the peptide or accession (i.e. protein) level. To retain protein IDs while aggregating to peptide level, set `aggregation_level <- c("accession","peptide")`. The entries are $log_2$-transformed after being normalized by the reference.


```r
# Set aggregation level
aggregation_level <- c("accession")
# Create cross-tab
crosstab <- create_crosstab(msnid, masic_data,
                            aggregation_level = aggregation_level,
                            fractions, samples, references)
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:20em; overflow-x: scroll; width:100%; "><table class="table" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">(\#tab:unnamed-chunk-21)<left>First 10 rows of the global quantitative cross-tab.</left>
</caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_126 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_127C </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_127N </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_128C </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_128N </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_129C </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_129N </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_130C </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_130N </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_126 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_127C </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_127N </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_128C </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_128N </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_129C </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_129N </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_130C </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_130N </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> AP_004893.1 </td>
   <td style="text-align:right;"> 0.1419768 </td>
   <td style="text-align:right;"> 0.7628195 </td>
   <td style="text-align:right;"> 0.1653552 </td>
   <td style="text-align:right;"> 0.8662554 </td>
   <td style="text-align:right;"> 0.9453172 </td>
   <td style="text-align:right;"> -0.6460065 </td>
   <td style="text-align:right;"> -1.9294467 </td>
   <td style="text-align:right;"> -0.4321433 </td>
   <td style="text-align:right;"> -1.2831873 </td>
   <td style="text-align:right;"> -1.0271227 </td>
   <td style="text-align:right;"> -0.9390945 </td>
   <td style="text-align:right;"> 0.4883309 </td>
   <td style="text-align:right;"> -1.7148628 </td>
   <td style="text-align:right;"> -0.7029685 </td>
   <td style="text-align:right;"> -0.8794712 </td>
   <td style="text-align:right;"> -0.1912097 </td>
   <td style="text-align:right;"> 0.3964607 </td>
   <td style="text-align:right;"> -0.2440478 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AP_004894.1 </td>
   <td style="text-align:right;"> 0.8092676 </td>
   <td style="text-align:right;"> -0.0976095 </td>
   <td style="text-align:right;"> -0.3113350 </td>
   <td style="text-align:right;"> 0.3215692 </td>
   <td style="text-align:right;"> 0.2171255 </td>
   <td style="text-align:right;"> -0.3678781 </td>
   <td style="text-align:right;"> -0.1638689 </td>
   <td style="text-align:right;"> -0.6696829 </td>
   <td style="text-align:right;"> -1.2039041 </td>
   <td style="text-align:right;"> -0.5124954 </td>
   <td style="text-align:right;"> -0.2364175 </td>
   <td style="text-align:right;"> -0.4428327 </td>
   <td style="text-align:right;"> -1.3730408 </td>
   <td style="text-align:right;"> -0.6711809 </td>
   <td style="text-align:right;"> -1.3515366 </td>
   <td style="text-align:right;"> -0.7462995 </td>
   <td style="text-align:right;"> -0.8338103 </td>
   <td style="text-align:right;"> -0.2227493 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AP_004895.1 </td>
   <td style="text-align:right;"> 0.2078433 </td>
   <td style="text-align:right;"> -0.2867209 </td>
   <td style="text-align:right;"> -0.6089756 </td>
   <td style="text-align:right;"> -0.1162062 </td>
   <td style="text-align:right;"> -0.3840271 </td>
   <td style="text-align:right;"> -1.1240967 </td>
   <td style="text-align:right;"> -0.6908468 </td>
   <td style="text-align:right;"> -0.6652575 </td>
   <td style="text-align:right;"> -0.7140383 </td>
   <td style="text-align:right;"> 0.2717217 </td>
   <td style="text-align:right;"> -0.1448289 </td>
   <td style="text-align:right;"> -0.1200736 </td>
   <td style="text-align:right;"> -0.6435709 </td>
   <td style="text-align:right;"> -0.4287771 </td>
   <td style="text-align:right;"> -0.6780284 </td>
   <td style="text-align:right;"> -0.6102404 </td>
   <td style="text-align:right;"> -0.3896190 </td>
   <td style="text-align:right;"> -0.1548544 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AP_004896.1 </td>
   <td style="text-align:right;"> -0.1494849 </td>
   <td style="text-align:right;"> -0.3664339 </td>
   <td style="text-align:right;"> -0.7314368 </td>
   <td style="text-align:right;"> -0.1742391 </td>
   <td style="text-align:right;"> -0.5352280 </td>
   <td style="text-align:right;"> -1.2945071 </td>
   <td style="text-align:right;"> -1.0372327 </td>
   <td style="text-align:right;"> -0.7060783 </td>
   <td style="text-align:right;"> -0.8299749 </td>
   <td style="text-align:right;"> 0.1939540 </td>
   <td style="text-align:right;"> -0.2274358 </td>
   <td style="text-align:right;"> -0.1688422 </td>
   <td style="text-align:right;"> -0.5251264 </td>
   <td style="text-align:right;"> -0.4222698 </td>
   <td style="text-align:right;"> -0.6543311 </td>
   <td style="text-align:right;"> -0.6741064 </td>
   <td style="text-align:right;"> -0.3994149 </td>
   <td style="text-align:right;"> -0.0441485 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AP_004898.1 </td>
   <td style="text-align:right;"> 0.0362964 </td>
   <td style="text-align:right;"> 0.7497227 </td>
   <td style="text-align:right;"> 0.4252227 </td>
   <td style="text-align:right;"> 0.4913660 </td>
   <td style="text-align:right;"> 1.1580326 </td>
   <td style="text-align:right;"> 0.1211536 </td>
   <td style="text-align:right;"> -0.3640632 </td>
   <td style="text-align:right;"> -0.3019505 </td>
   <td style="text-align:right;"> -0.8291744 </td>
   <td style="text-align:right;"> -0.8407749 </td>
   <td style="text-align:right;"> -0.2796091 </td>
   <td style="text-align:right;"> -0.4130732 </td>
   <td style="text-align:right;"> -1.5747761 </td>
   <td style="text-align:right;"> -0.9449498 </td>
   <td style="text-align:right;"> -1.8439756 </td>
   <td style="text-align:right;"> -0.1774225 </td>
   <td style="text-align:right;"> -1.1083199 </td>
   <td style="text-align:right;"> -0.4175363 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AP_004899.1 </td>
   <td style="text-align:right;"> 0.7140968 </td>
   <td style="text-align:right;"> -0.1781542 </td>
   <td style="text-align:right;"> -0.3732752 </td>
   <td style="text-align:right;"> 0.3494902 </td>
   <td style="text-align:right;"> -0.0615626 </td>
   <td style="text-align:right;"> -2.1679002 </td>
   <td style="text-align:right;"> -0.8550940 </td>
   <td style="text-align:right;"> -0.9026145 </td>
   <td style="text-align:right;"> -1.4519278 </td>
   <td style="text-align:right;"> -0.3158081 </td>
   <td style="text-align:right;"> -0.4056811 </td>
   <td style="text-align:right;"> -0.4644758 </td>
   <td style="text-align:right;"> -0.2805080 </td>
   <td style="text-align:right;"> -0.9023044 </td>
   <td style="text-align:right;"> -1.0482424 </td>
   <td style="text-align:right;"> -0.8052899 </td>
   <td style="text-align:right;"> -0.6675429 </td>
   <td style="text-align:right;"> -0.3959923 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AP_004900.1 </td>
   <td style="text-align:right;"> -0.3806966 </td>
   <td style="text-align:right;"> -0.3441177 </td>
   <td style="text-align:right;"> -0.5883203 </td>
   <td style="text-align:right;"> -0.0902205 </td>
   <td style="text-align:right;"> -0.8263700 </td>
   <td style="text-align:right;"> -0.7060111 </td>
   <td style="text-align:right;"> -1.0978191 </td>
   <td style="text-align:right;"> -0.8570849 </td>
   <td style="text-align:right;"> -1.0769673 </td>
   <td style="text-align:right;"> -0.1566909 </td>
   <td style="text-align:right;"> -0.2565750 </td>
   <td style="text-align:right;"> -0.5707603 </td>
   <td style="text-align:right;"> -0.5960161 </td>
   <td style="text-align:right;"> -0.6380722 </td>
   <td style="text-align:right;"> -0.5524057 </td>
   <td style="text-align:right;"> -0.6422737 </td>
   <td style="text-align:right;"> -0.5140577 </td>
   <td style="text-align:right;"> -0.1386396 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AP_004901.1 </td>
   <td style="text-align:right;"> -0.2839471 </td>
   <td style="text-align:right;"> -0.5758177 </td>
   <td style="text-align:right;"> -0.0216835 </td>
   <td style="text-align:right;"> -0.2468966 </td>
   <td style="text-align:right;"> 0.3802436 </td>
   <td style="text-align:right;"> -0.7852805 </td>
   <td style="text-align:right;"> -1.2036962 </td>
   <td style="text-align:right;"> -1.0623455 </td>
   <td style="text-align:right;"> -1.3738888 </td>
   <td style="text-align:right;"> 0.0884387 </td>
   <td style="text-align:right;"> -0.5271134 </td>
   <td style="text-align:right;"> -0.3024421 </td>
   <td style="text-align:right;"> -0.8704195 </td>
   <td style="text-align:right;"> -0.5130666 </td>
   <td style="text-align:right;"> -0.8245858 </td>
   <td style="text-align:right;"> -0.3465590 </td>
   <td style="text-align:right;"> -0.3041638 </td>
   <td style="text-align:right;"> -0.0778906 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NP_001000613.1 </td>
   <td style="text-align:right;"> 0.0642942 </td>
   <td style="text-align:right;"> -0.2044263 </td>
   <td style="text-align:right;"> -0.4935404 </td>
   <td style="text-align:right;"> 0.1065480 </td>
   <td style="text-align:right;"> -0.3240158 </td>
   <td style="text-align:right;"> 0.0233030 </td>
   <td style="text-align:right;"> -0.7883731 </td>
   <td style="text-align:right;"> -0.5831656 </td>
   <td style="text-align:right;"> -1.1649507 </td>
   <td style="text-align:right;"> -0.0969739 </td>
   <td style="text-align:right;"> -0.2298260 </td>
   <td style="text-align:right;"> -0.5051282 </td>
   <td style="text-align:right;"> -0.9489927 </td>
   <td style="text-align:right;"> -0.7860575 </td>
   <td style="text-align:right;"> -0.4803511 </td>
   <td style="text-align:right;"> -0.4824157 </td>
   <td style="text-align:right;"> -0.8430505 </td>
   <td style="text-align:right;"> -0.4220978 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NP_001001512.2 </td>
   <td style="text-align:right;"> -0.2155278 </td>
   <td style="text-align:right;"> -0.1288367 </td>
   <td style="text-align:right;"> -0.4900798 </td>
   <td style="text-align:right;"> -0.2239173 </td>
   <td style="text-align:right;"> -0.4303687 </td>
   <td style="text-align:right;"> -0.2496957 </td>
   <td style="text-align:right;"> -0.9165473 </td>
   <td style="text-align:right;"> -0.4952781 </td>
   <td style="text-align:right;"> -1.0160964 </td>
   <td style="text-align:right;"> -0.1256882 </td>
   <td style="text-align:right;"> -0.2742136 </td>
   <td style="text-align:right;"> -0.2713873 </td>
   <td style="text-align:right;"> -0.5914317 </td>
   <td style="text-align:right;"> -0.4021118 </td>
   <td style="text-align:right;"> -0.8604247 </td>
   <td style="text-align:right;"> -0.6648978 </td>
   <td style="text-align:right;"> -0.6082902 </td>
   <td style="text-align:right;"> -0.1742713 </td>
  </tr>
</tbody>
</table></div>

</br>

In order to demonstrate prioritized inference in Section \@ref(phosphoproteomics-data), we need to save the row names of this cross-tab.


```r
# Save protein names
saveRDS(rownames(crosstab), file = "data/3442_global_protein_names.rds")
```

We should save the cross-tab as well. To do so, we need to convert the row names to a column called `protein`.


```r
# Modify cross-tab for saving
crosstab <- crosstab %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("protein")

# Save cross-tab
write.table(crosstab, file = "data/global_quant_crosstab.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
```


## Phosphoproteomics Data {#phosphoproteomics-data}

This pipeline shows how to process data from the DMS. The number of the data package is `3626`. For this section, we need the `PlexedPiper` package for isobaric quantification and `PNNL.DMS.utils` to interface with the DMS. Also, some details will be omitted if they were already provided in Section \@ref(global-proteomics-data).


```r
# Setup
library(PNNL.DMS.utils)
library(PlexedPiper)
library(Biostrings)
library(dplyr) # %>%
```


### Prepare MS/MS Identifications {#prepare-MS2-IDs-phospho}

#### Read MS-GF+ Output {-}


```r
# Read MS-GF+ data
data_package_num <- 3626
msnid <- read_msgf_data_from_DMS(data_package_num)
```


```r
show(msnid)
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  23 
## #PSMs: 612667 at 55 % FDR
## #peptides: 396540 at 75 % FDR
## #accessions: 121521 at 98 % FDR
```


#### Remove Non-Phosphorylated Peptides {-}

In this case, the phosphorylation of an amino acid is marked by a `*` appearing after the amino acid. We will not consider unmodified peptides, so we can filter them out. The `*` is a special character that must be escaped with backslashes, and the backslashes must also be escaped.


```r
# Remove non-phosphorylated peptides
# (peptides that do not contain a *)
msnid <- apply_filter(msnid, "grepl('\\\\*', peptide)")
show(msnid)
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  23 
## #PSMs: 537749 at 57 % FDR
## #peptides: 353634 at 76 % FDR
## #accessions: 118817 at 98 % FDR
```


#### Correct Isotope Selection Error {-}


```r
# Correct for isotope selection error
msnid <- correct_peak_selection(msnid)
```


#### Remove Contaminants {-}


```r
# Remove contaminants
msnid <- apply_filter(msnid, "!grepl('Contaminant', accession)")
show(msnid)
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  23 
## #PSMs: 537572 at 57 % FDR
## #peptides: 353489 at 76 % FDR
## #accessions: 118797 at 98 % FDR
```


#### AScore {-}

Phospho datasets involve AScore jobs for improving phosphosite localization. There should be one AScore job per data package. The fetched object is a data.frame that links datasets, scans and original PTM localization to newly suggested locations. Importantly, it contains `AScore` column that signifies the confidence of PTM assignment. AScore > 17 is considered confident.


```r
# Filter PTMs by AScore
ascore <- get_AScore_results(data_package_num)
msnid <- best_PTM_location_by_ascore(msnid, ascore)
```


```r
show(msnid)
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  23 
## #PSMs: 188791 at 30 % FDR
## #peptides: 101873 at 53 % FDR
## #accessions: 90677 at 93 % FDR
```


#### MS/MS ID Filter: Peptide Level {-}


```r
# 1% FDR filter at the peptide level
msnid <- filter_msgf_data(msnid,
                          level = "peptide",
                          fdr.max = 0.01)
show(msnid)
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  23 
## #PSMs: 76103 at 0.49 % FDR
## #peptides: 23378 at 1 % FDR
## #accessions: 16090 at 4.7 % FDR
```

#### MS/MS ID Filter: Protein Level {-}


```r
# Get path to FASTA file
path_to_FASTA <- path_to_FASTA_used_by_DMS(data_package_num)

# Compute number of peptides per 1000 amino acids
msnid <- compute_num_peptides_per_1000aa(msnid, path_to_FASTA)

# 1% FDR filter at the protein level
msnid <- filter_msgf_data(msnid,
                          level = "accession",
                          fdr.max = 0.01)
show(msnid)
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  23 
## #PSMs: 72471 at 0.12 % FDR
## #peptides: 21261 at 0.25 % FDR
## #accessions: 9413 at 0.95 % FDR
```

#### Inference of Parsimonious Protein Set {-}

<!---
TODO:
* Talk about using prior information from global cross-tab to improve inference of the parsimonious protein set for phospho data.
--->


```r
# Load proteins from global crosstab
global_proteins <- readRDS("data/3442_global_protein_names.rds")
# Inference of parsimonious protein set
msnid <- infer_parsimonious_accessions(msnid, unique_only = FALSE,
                                       prior = global_proteins)
show(msnid)
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  23 
## #PSMs: 72471 at 0.12 % FDR
## #peptides: 21261 at 0.25 % FDR
## #accessions: 2890 at 1.6 % FDR
```

Notice that the protein-level FDR increased above the acceptable threshold, so we need to reapply the filter.


```r
# 1% FDR filter at the protein level
msnid <- filter_msgf_data(msnid,
                          level = "accession",
                          fdr.max = 0.01)
show(msnid)
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  23 
## #PSMs: 70359 at 0.075 % FDR
## #peptides: 20127 at 0.15 % FDR
## #accessions: 2356 at 0.99 % FDR
```


```r
# Remove Decoy Accessions
msnid <- apply_filter(msnid, "!isDecoy")
show(msnid)
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  23 
## #PSMs: 70306 at 0 % FDR
## #peptides: 20096 at 0 % FDR
## #accessions: 2333 at 0 % FDR
```

#### Map Sites to Protein Sequences {-}

Prepare FASTA to make sure entry names in FASTA file match MSnID accessions. The plan is to make this conversion automatic. `map_mod_sites` creates number of columns describing mapping of the site/s onto the protein sequences. The most important for the user is `SiteID`.


```r
# Create AAStringSet
fst <- readAAStringSet(path_to_FASTA)
# Remove contaminants
fst <- fst[!grepl("Contaminant", names(fst)), ]
# First 6 names
head(names(fst))
```

```
## [1] "NP_783171.2 cathepsin R precursor [Rattus norvegicus]"                    
## [2] "NP_001101862.2 zinc finger protein ZIC 2 [Rattus norvegicus]"             
## [3] "NP_113721.4 UDP-glucuronosyltransferase 2B2 precursor [Rattus norvegicus]"
## [4] "NP_714948.1 Ly-49 stimulatory receptor 3 [Rattus norvegicus]"             
## [5] "NP_001000704.1 olfactory receptor Olr931 [Rattus norvegicus]"             
## [6] "NP_001000638.1 olfactory receptor Olr652 [Rattus norvegicus]"
```


```r
# Modify names to match accessions(msnid)
names(fst) <- strsplit(names(fst), split = " ") %>% 
  # Select text before first space
  lapply(function(x) x[1]) %>% 
  unlist()
# First 6 names
head(names(fst))
```

```
## [1] "NP_783171.2"    "NP_001101862.2" "NP_113721.4"    "NP_714948.1"   
## [5] "NP_001000704.1" "NP_001000638.1"
```


```r
# Main mapping call
msnid <- map_mod_sites(object = msnid, fasta = fst, 
                       accession_col = "accession", 
                       peptide_mod_col = "peptide", 
                       mod_char = "*",
                       site_delimiter = "lower")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:20em; overflow-x: scroll; width:100%; "><table class="table" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">(\#tab:phospho-msnid-table)<left>First 10 rows of the processed MS-GF+ results.</left>
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
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> OriginalPeptide </th>
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
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> peptide </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> maxAScore </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> msmsScore </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> absParentMassErrorPPM </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> peptides_per_1000aa </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> First_AA </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Last_AA </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> First_AA_First </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> Last_AA_First </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> ProtLen </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> ModShift </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> ModAAs </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SiteLoc </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> Site </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SiteCollapsed </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SiteCollapsedFirst </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;"> SiteID </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_P_S1_06_DIL_28Oct17_Elm_AQ-17-10-03 </td>
   <td style="text-align:right;"> 12697 </td>
   <td style="text-align:right;"> 27321 </td>
   <td style="text-align:left;"> HCD </td>
   <td style="text-align:right;"> 2256 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 1045.124 </td>
   <td style="text-align:right;"> 0.003 </td>
   <td style="text-align:right;"> 0.858 </td>
   <td style="text-align:right;"> 3131.346 </td>
   <td style="text-align:left;"> A.AAAAAGDS\*DS\*WDADTFSMEDPVRK.V </td>
   <td style="text-align:left;"> NP_001071138.1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 146 </td>
   <td style="text-align:right;"> 58 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> NP_001071138.1 </td>
   <td style="text-align:right;"> 1044.454 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 1044.455 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_P_S1_06_DIL_28Oct17_Elm_AQ-17-10-03 </td>
   <td style="text-align:right;"> 27321 </td>
   <td style="text-align:left;"> AAAAAGDSDSWDADTFSMEDPVRK </td>
   <td style="text-align:left;"> A.AAAAAGDSDS\*WDADT\*FSMEDPVRK.V </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> 1.057 </td>
   <td style="text-align:right;"> 7.722 </td>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:left;"> 28 </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 28 </td>
   <td style="text-align:right;"> 259 </td>
   <td style="text-align:left;"> 9, 14 </td>
   <td style="text-align:left;"> S, T </td>
   <td style="text-align:left;"> 14, 19 </td>
   <td style="text-align:left;"> S14, T19 </td>
   <td style="text-align:left;"> S14,T19 </td>
   <td style="text-align:left;"> S14,T19 </td>
   <td style="text-align:left;"> NP_001071138.1-S14sT19t </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_P_S1_07_DIL_28Oct17_Elm_AQ-17-10-03 </td>
   <td style="text-align:right;"> 875 </td>
   <td style="text-align:right;"> 23519 </td>
   <td style="text-align:left;"> HCD </td>
   <td style="text-align:right;"> 264 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 952.144 </td>
   <td style="text-align:right;"> 0.004 </td>
   <td style="text-align:right;"> 1.538 </td>
   <td style="text-align:right;"> 2854.412 </td>
   <td style="text-align:left;"> R.AAAASAAEAGIAT\*PGTEDSDDALLK.M </td>
   <td style="text-align:left;"> XP_006232986.1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 165 </td>
   <td style="text-align:right;"> 129 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> XP_006232986.1 </td>
   <td style="text-align:right;"> 952.142 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 952.144 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_P_S1_07_DIL_28Oct17_Elm_AQ-17-10-03 </td>
   <td style="text-align:right;"> 23519 </td>
   <td style="text-align:left;"> AAAASAAEAGIATPGTEDSDDALLK </td>
   <td style="text-align:left;"> R.AAAASAAEAGIAT\*PGTEDSDDALLK.M </td>
   <td style="text-align:right;"> 52.349 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> 1.625 </td>
   <td style="text-align:right;"> 5.305 </td>
   <td style="text-align:left;"> 238 </td>
   <td style="text-align:left;"> 262 </td>
   <td style="text-align:right;"> 238 </td>
   <td style="text-align:right;"> 262 </td>
   <td style="text-align:right;"> 377 </td>
   <td style="text-align:left;"> 12 </td>
   <td style="text-align:left;"> T </td>
   <td style="text-align:left;"> 250 </td>
   <td style="text-align:left;"> T250 </td>
   <td style="text-align:left;"> T250 </td>
   <td style="text-align:left;"> T250 </td>
   <td style="text-align:left;"> XP_006232986.1-T250t </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_P_S1_07_DIL_28Oct17_Elm_AQ-17-10-03 </td>
   <td style="text-align:right;"> 12873 </td>
   <td style="text-align:right;"> 23508 </td>
   <td style="text-align:left;"> HCD </td>
   <td style="text-align:right;"> 2213 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 714.360 </td>
   <td style="text-align:right;"> 0.007 </td>
   <td style="text-align:right;"> 2.392 </td>
   <td style="text-align:right;"> 2854.412 </td>
   <td style="text-align:left;"> R.AAAASAAEAGIAT\*PGTEDSDDALLK.M </td>
   <td style="text-align:left;"> XP_006232986.1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 122 </td>
   <td style="text-align:right;"> 81 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:left;"> XP_006232986.1 </td>
   <td style="text-align:right;"> 714.358 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 714.360 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_P_S1_07_DIL_28Oct17_Elm_AQ-17-10-03 </td>
   <td style="text-align:right;"> 23508 </td>
   <td style="text-align:left;"> AAAASAAEAGIATPGTEDSDDALLK </td>
   <td style="text-align:left;"> R.AAAASAAEAGIAT\*PGTEDSDDALLK.M </td>
   <td style="text-align:right;"> 17.480 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> 2.472 </td>
   <td style="text-align:right;"> 5.305 </td>
   <td style="text-align:left;"> 238 </td>
   <td style="text-align:left;"> 262 </td>
   <td style="text-align:right;"> 238 </td>
   <td style="text-align:right;"> 262 </td>
   <td style="text-align:right;"> 377 </td>
   <td style="text-align:left;"> 12 </td>
   <td style="text-align:left;"> T </td>
   <td style="text-align:left;"> 250 </td>
   <td style="text-align:left;"> T250 </td>
   <td style="text-align:left;"> T250 </td>
   <td style="text-align:left;"> T250 </td>
   <td style="text-align:left;"> XP_006232986.1-T250t </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_P_S2_07_3Nov17_Elm_AQ-17-10-03 </td>
   <td style="text-align:right;"> 2731 </td>
   <td style="text-align:right;"> 23697 </td>
   <td style="text-align:left;"> HCD </td>
   <td style="text-align:right;"> 502 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 714.610 </td>
   <td style="text-align:right;"> 0.002 </td>
   <td style="text-align:right;"> 0.706 </td>
   <td style="text-align:right;"> 2854.412 </td>
   <td style="text-align:left;"> R.AAAASAAEAGIAT\*PGTEDSDDALLK.M </td>
   <td style="text-align:left;"> XP_006232986.1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 135 </td>
   <td style="text-align:right;"> 104 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> XP_006232986.1 </td>
   <td style="text-align:right;"> 714.358 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 714.359 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_P_S2_07_3Nov17_Elm_AQ-17-10-03 </td>
   <td style="text-align:right;"> 23697 </td>
   <td style="text-align:left;"> AAAASAAEAGIATPGTEDSDDALLK </td>
   <td style="text-align:left;"> R.AAAASAAEAGIAT\*PGTEDSDDALLK.M </td>
   <td style="text-align:right;"> 26.295 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> 0.780 </td>
   <td style="text-align:right;"> 5.305 </td>
   <td style="text-align:left;"> 238 </td>
   <td style="text-align:left;"> 262 </td>
   <td style="text-align:right;"> 238 </td>
   <td style="text-align:right;"> 262 </td>
   <td style="text-align:right;"> 377 </td>
   <td style="text-align:left;"> 12 </td>
   <td style="text-align:left;"> T </td>
   <td style="text-align:left;"> 250 </td>
   <td style="text-align:left;"> T250 </td>
   <td style="text-align:left;"> T250 </td>
   <td style="text-align:left;"> T250 </td>
   <td style="text-align:left;"> XP_006232986.1-T250t </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_P_S1_07_DIL_28Oct17_Elm_AQ-17-10-03 </td>
   <td style="text-align:right;"> 4877 </td>
   <td style="text-align:right;"> 21265 </td>
   <td style="text-align:left;"> HCD </td>
   <td style="text-align:right;"> 935 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 800.403 </td>
   <td style="text-align:right;"> 0.006 </td>
   <td style="text-align:right;"> 1.871 </td>
   <td style="text-align:right;"> 3196.577 </td>
   <td style="text-align:left;"> R.AAAASAAEAGIAT\*PGTEGERDSDDALLK.M </td>
   <td style="text-align:left;"> NP_112621.1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 194 </td>
   <td style="text-align:right;"> 114 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> NP_112621.1 </td>
   <td style="text-align:right;"> 799.900 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 799.901 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_P_S1_07_DIL_28Oct17_Elm_AQ-17-10-03 </td>
   <td style="text-align:right;"> 21265 </td>
   <td style="text-align:left;"> AAAASAAEAGIATPGTEGERDSDDALLK </td>
   <td style="text-align:left;"> R.AAAASAAEAGIATPGT\*EGERDSDDALLK.M </td>
   <td style="text-align:right;"> 6.213 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> 1.902 </td>
   <td style="text-align:right;"> 10.526 </td>
   <td style="text-align:left;"> 238 </td>
   <td style="text-align:left;"> 265 </td>
   <td style="text-align:right;"> 238 </td>
   <td style="text-align:right;"> 265 </td>
   <td style="text-align:right;"> 380 </td>
   <td style="text-align:left;"> 15 </td>
   <td style="text-align:left;"> T </td>
   <td style="text-align:left;"> 253 </td>
   <td style="text-align:left;"> T253 </td>
   <td style="text-align:left;"> T253 </td>
   <td style="text-align:left;"> T253 </td>
   <td style="text-align:left;"> NP_112621.1-T253t </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_P_S1_07_DIL_28Oct17_Elm_AQ-17-10-03 </td>
   <td style="text-align:right;"> 6826 </td>
   <td style="text-align:right;"> 21280 </td>
   <td style="text-align:left;"> HCD </td>
   <td style="text-align:right;"> 1251 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 1066.532 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> -0.095 </td>
   <td style="text-align:right;"> 3196.577 </td>
   <td style="text-align:left;"> R.AAAASAAEAGIATPGT\*EGERDSDDALLK.M </td>
   <td style="text-align:left;"> NP_112621.1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 200 </td>
   <td style="text-align:right;"> 94 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NP_112621.1 </td>
   <td style="text-align:right;"> 1066.197 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 1066.197 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_P_S1_07_DIL_28Oct17_Elm_AQ-17-10-03 </td>
   <td style="text-align:right;"> 21280 </td>
   <td style="text-align:left;"> AAAASAAEAGIATPGTEGERDSDDALLK </td>
   <td style="text-align:left;"> R.AAAASAAEAGIATPGT\*EGERDSDDALLK.M </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> 0.043 </td>
   <td style="text-align:right;"> 10.526 </td>
   <td style="text-align:left;"> 238 </td>
   <td style="text-align:left;"> 265 </td>
   <td style="text-align:right;"> 238 </td>
   <td style="text-align:right;"> 265 </td>
   <td style="text-align:right;"> 380 </td>
   <td style="text-align:left;"> 15 </td>
   <td style="text-align:left;"> T </td>
   <td style="text-align:left;"> 253 </td>
   <td style="text-align:left;"> T253 </td>
   <td style="text-align:left;"> T253 </td>
   <td style="text-align:left;"> T253 </td>
   <td style="text-align:left;"> NP_112621.1-T253t </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_P_S1_09_DIL_28Oct17_Elm_AQ-17-10-03 </td>
   <td style="text-align:right;"> 4625 </td>
   <td style="text-align:right;"> 20810 </td>
   <td style="text-align:left;"> HCD </td>
   <td style="text-align:right;"> 736 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 800.150 </td>
   <td style="text-align:right;"> -0.002 </td>
   <td style="text-align:right;"> -0.743 </td>
   <td style="text-align:right;"> 3196.577 </td>
   <td style="text-align:left;"> R.AAAASAAEAGIATPGT\*EGERDSDDALLK.M </td>
   <td style="text-align:left;"> NP_112621.1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 156 </td>
   <td style="text-align:right;"> 93 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NP_112621.1 </td>
   <td style="text-align:right;"> 799.900 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;"> 799.899 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_P_S1_09_DIL_28Oct17_Elm_AQ-17-10-03 </td>
   <td style="text-align:right;"> 20810 </td>
   <td style="text-align:left;"> AAAASAAEAGIATPGTEGERDSDDALLK </td>
   <td style="text-align:left;"> R.AAAASAAEAGIATPGTEGERDS\*DDALLK.M </td>
   <td style="text-align:right;"> 5.771 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> 0.712 </td>
   <td style="text-align:right;"> 10.526 </td>
   <td style="text-align:left;"> 238 </td>
   <td style="text-align:left;"> 265 </td>
   <td style="text-align:right;"> 238 </td>
   <td style="text-align:right;"> 265 </td>
   <td style="text-align:right;"> 380 </td>
   <td style="text-align:left;"> 21 </td>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> 259 </td>
   <td style="text-align:left;"> S259 </td>
   <td style="text-align:left;"> S259 </td>
   <td style="text-align:left;"> S259 </td>
   <td style="text-align:left;"> NP_112621.1-S259s </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_P_S1_09_DIL_28Oct17_Elm_AQ-17-10-03 </td>
   <td style="text-align:right;"> 16605 </td>
   <td style="text-align:right;"> 20839 </td>
   <td style="text-align:left;"> HCD </td>
   <td style="text-align:right;"> 2489 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 1066.530 </td>
   <td style="text-align:right;"> -0.006 </td>
   <td style="text-align:right;"> -1.812 </td>
   <td style="text-align:right;"> 3196.577 </td>
   <td style="text-align:left;"> R.AAAASAAEAGIATPGTEGERDS\*DDALLK.M </td>
   <td style="text-align:left;"> NP_112621.1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 159 </td>
   <td style="text-align:right;"> 48 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.01 </td>
   <td style="text-align:right;"> 0.005 </td>
   <td style="text-align:right;"> 0.007 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NP_112621.1 </td>
   <td style="text-align:right;"> 1066.197 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 1066.195 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_P_S1_09_DIL_28Oct17_Elm_AQ-17-10-03 </td>
   <td style="text-align:right;"> 20839 </td>
   <td style="text-align:left;"> AAAASAAEAGIATPGTEGERDSDDALLK </td>
   <td style="text-align:left;"> R.AAAASAAEAGIATPGTEGERDS\*DDALLK.M </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 2.148 </td>
   <td style="text-align:right;"> 1.673 </td>
   <td style="text-align:right;"> 10.526 </td>
   <td style="text-align:left;"> 238 </td>
   <td style="text-align:left;"> 265 </td>
   <td style="text-align:right;"> 238 </td>
   <td style="text-align:right;"> 265 </td>
   <td style="text-align:right;"> 380 </td>
   <td style="text-align:left;"> 21 </td>
   <td style="text-align:left;"> S </td>
   <td style="text-align:left;"> 259 </td>
   <td style="text-align:left;"> S259 </td>
   <td style="text-align:left;"> S259 </td>
   <td style="text-align:left;"> S259 </td>
   <td style="text-align:left;"> NP_112621.1-S259s </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_P_S2_07_3Nov17_Elm_AQ-17-10-03 </td>
   <td style="text-align:right;"> 433 </td>
   <td style="text-align:right;"> 21424 </td>
   <td style="text-align:left;"> HCD </td>
   <td style="text-align:right;"> 54 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 1066.531 </td>
   <td style="text-align:right;"> -0.004 </td>
   <td style="text-align:right;"> -1.125 </td>
   <td style="text-align:right;"> 3196.577 </td>
   <td style="text-align:left;"> R.AAAASAAEAGIAT\*PGTEGERDSDDALLK.M </td>
   <td style="text-align:left;"> NP_112621.1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 241 </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NP_112621.1 </td>
   <td style="text-align:right;"> 1066.197 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 1066.196 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_P_S2_07_3Nov17_Elm_AQ-17-10-03 </td>
   <td style="text-align:right;"> 21424 </td>
   <td style="text-align:left;"> AAAASAAEAGIATPGTEGERDSDDALLK </td>
   <td style="text-align:left;"> R.AAAASAAEAGIAT\*PGTEGERDSDDALLK.M </td>
   <td style="text-align:right;"> 32.347 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> 0.989 </td>
   <td style="text-align:right;"> 10.526 </td>
   <td style="text-align:left;"> 238 </td>
   <td style="text-align:left;"> 265 </td>
   <td style="text-align:right;"> 238 </td>
   <td style="text-align:right;"> 265 </td>
   <td style="text-align:right;"> 380 </td>
   <td style="text-align:left;"> 12 </td>
   <td style="text-align:left;"> T </td>
   <td style="text-align:left;"> 250 </td>
   <td style="text-align:left;"> T250 </td>
   <td style="text-align:left;"> T250 </td>
   <td style="text-align:left;"> T250 </td>
   <td style="text-align:left;"> NP_112621.1-T250t </td>
  </tr>
  <tr>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_P_S2_07_3Nov17_Elm_AQ-17-10-03 </td>
   <td style="text-align:right;"> 434 </td>
   <td style="text-align:right;"> 21424 </td>
   <td style="text-align:left;"> HCD </td>
   <td style="text-align:right;"> 54 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 1066.531 </td>
   <td style="text-align:right;"> -0.004 </td>
   <td style="text-align:right;"> -1.125 </td>
   <td style="text-align:right;"> 3196.577 </td>
   <td style="text-align:left;"> R.AAAASAAEAGIATPGT\*EGERDSDDALLK.M </td>
   <td style="text-align:left;"> NP_112621.1 </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 241 </td>
   <td style="text-align:right;"> 168 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 0.000 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> NP_112621.1 </td>
   <td style="text-align:right;"> 1066.197 </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 1066.196 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> MoTrPAC_Pilot_TMT_P_S2_07_3Nov17_Elm_AQ-17-10-03 </td>
   <td style="text-align:right;"> 21424 </td>
   <td style="text-align:left;"> AAAASAAEAGIATPGTEGERDSDDALLK </td>
   <td style="text-align:left;"> R.AAAASAAEAGIAT\*PGTEGERDSDDALLK.M </td>
   <td style="text-align:right;"> 32.347 </td>
   <td style="text-align:right;"> Inf </td>
   <td style="text-align:right;"> 0.989 </td>
   <td style="text-align:right;"> 10.526 </td>
   <td style="text-align:left;"> 238 </td>
   <td style="text-align:left;"> 265 </td>
   <td style="text-align:right;"> 238 </td>
   <td style="text-align:right;"> 265 </td>
   <td style="text-align:right;"> 380 </td>
   <td style="text-align:left;"> 12 </td>
   <td style="text-align:left;"> T </td>
   <td style="text-align:left;"> 250 </td>
   <td style="text-align:left;"> T250 </td>
   <td style="text-align:left;"> T250 </td>
   <td style="text-align:left;"> T250 </td>
   <td style="text-align:left;"> NP_112621.1-T250t </td>
  </tr>
</tbody>
</table></div>

</br>

Table \@ref(tab:phospho-msnid-table) shows the first 10 rows of the processed MS-GF+ output.


### Prepare Reporter Ion Intensities {#phospho-reporter-ion-intensities}

#### Read MASIC Output {-}


```r
# Read MASIC data
masic_data <- read_masic_data_from_DMS(data_package_num, 
                                       interference_score = TRUE)
```


#### Filter MASIC Data {-}


```r
# Filter MASIC data
masic_data <- filter_masic_data(masic_data, 
                                interference_score_threshold = 0.5,
                                s2n_threshold = 0)
```


### Create Study Design Tables

#### Fractions {-}


```r
# Create fractions table
fractions <- data.frame(Dataset = unique(masic_data$Dataset)) %>% 
  mutate(PlexID = gsub(".*_P_(S\\d{1})_.*", "\\1", Dataset))
```


#### Samples {-}


```r
# TMT10 Reporter Converter table from MSnID package
conv <- reporter_converter$tmt10
plexes <- unique(fractions$PlexID)

# Reference channel
ref_channel <- "131"

# Create samples table
samples <- data.frame(PlexID = rep(plexes, each = nrow(conv)),  
                      ReporterName = rep(conv$ReporterName, 
                                         length(plexes))) %>% 
  mutate(ReporterAlias = sprintf("%s_%s", PlexID, ReporterName),
         MeasurementName = ReporterAlias,
         QuantBlock = 1,
         # Comment out this next part if the reference
         # is not one of the reporter ion channels.
         MeasurementName = ifelse(ReporterName == ref_channel,
                                  NA, MeasurementName)
  )
```

#### References {-}


```r
# Create references table
references <- samples %>% 
  # Filter to reference channel
  filter(ReporterName == ref_channel) %>% 
  # Select required columns and rename ReporterAlias to Reference
  select(PlexID, Reference = ReporterAlias, QuantBlock)
```


### Create Quantitative Cross-tab {#phospho-create-crosstab}


```r
# Set aggregation level
aggregation_level <- c("accession", "peptide")
# Create cross-tab
crosstab <- create_crosstab(msnid, masic_data,
                            aggregation_level = aggregation_level,
                            fractions, samples, references)
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:20em; overflow-x: scroll; width:100%; "><table class="table" style="font-size: 12px; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">(\#tab:unnamed-chunk-46)<left>First 10 rows of the phospho quantitative cross-tab.</left>
</caption>
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_126 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_127C </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_127N </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_128C </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_128N </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_129C </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_129N </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_130C </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S1_130N </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_126 </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_127C </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_127N </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_128C </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_128N </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_129C </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_129N </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_130C </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;"> S2_130N </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> NP_001001512.2@\hphantom{}K.S\*SEPPPPPPVPEPTNAGK.R </td>
   <td style="text-align:right;"> -0.5441083 </td>
   <td style="text-align:right;"> -0.0291611 </td>
   <td style="text-align:right;"> -0.5540885 </td>
   <td style="text-align:right;"> -0.1812456 </td>
   <td style="text-align:right;"> -0.5675833 </td>
   <td style="text-align:right;"> -0.0427948 </td>
   <td style="text-align:right;"> -0.8030560 </td>
   <td style="text-align:right;"> -0.5236890 </td>
   <td style="text-align:right;"> -1.0398773 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NP_001001512.2@\hphantom{}K.SS\*EPPPPPPVPEPTNAGK.R </td>
   <td style="text-align:right;"> -0.2806018 </td>
   <td style="text-align:right;"> -0.2939827 </td>
   <td style="text-align:right;"> -0.3923442 </td>
   <td style="text-align:right;"> -0.3242804 </td>
   <td style="text-align:right;"> -0.6512914 </td>
   <td style="text-align:right;"> -0.4914234 </td>
   <td style="text-align:right;"> -1.2105927 </td>
   <td style="text-align:right;"> -0.6421375 </td>
   <td style="text-align:right;"> -0.8869142 </td>
   <td style="text-align:right;"> -0.3079991 </td>
   <td style="text-align:right;"> -0.8072458 </td>
   <td style="text-align:right;"> -0.4801256 </td>
   <td style="text-align:right;"> -0.4458455 </td>
   <td style="text-align:right;"> -0.7588038 </td>
   <td style="text-align:right;"> -1.3958746 </td>
   <td style="text-align:right;"> -0.6640721 </td>
   <td style="text-align:right;"> -0.7046127 </td>
   <td style="text-align:right;"> -0.3091972 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NP_001001512.2@\hphantom{}K.TNSS\*PSVNTTASGVEDLNIIQVTIPDDDNER.L </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -1.4611181 </td>
   <td style="text-align:right;"> -2.2975481 </td>
   <td style="text-align:right;"> -1.6828088 </td>
   <td style="text-align:right;"> -2.3039510 </td>
   <td style="text-align:right;"> -1.8632844 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -1.5981389 </td>
   <td style="text-align:right;"> -2.4732843 </td>
   <td style="text-align:right;"> -0.4736739 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NP_001001512.2@\hphantom{}K.TNSSPS\*VNTTASGVEDLNIIQVTIPDDDNER.L </td>
   <td style="text-align:right;"> -1.2580276 </td>
   <td style="text-align:right;"> -1.3182692 </td>
   <td style="text-align:right;"> -2.1938861 </td>
   <td style="text-align:right;"> -2.0017857 </td>
   <td style="text-align:right;"> 0.0219220 </td>
   <td style="text-align:right;"> -0.7086501 </td>
   <td style="text-align:right;"> -1.6940153 </td>
   <td style="text-align:right;"> -1.4681967 </td>
   <td style="text-align:right;"> -1.5795032 </td>
   <td style="text-align:right;"> 0.3139811 </td>
   <td style="text-align:right;"> -0.6805752 </td>
   <td style="text-align:right;"> -0.8132341 </td>
   <td style="text-align:right;"> -1.7933771 </td>
   <td style="text-align:right;"> 0.0127816 </td>
   <td style="text-align:right;"> 0.2726023 </td>
   <td style="text-align:right;"> 0.2266855 </td>
   <td style="text-align:right;"> 0.5058189 </td>
   <td style="text-align:right;"> 0.9332038 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NP_001001512.2@\hphantom{}R.RPS\*TFGIPR.L </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -0.7930334 </td>
   <td style="text-align:right;"> -1.2173509 </td>
   <td style="text-align:right;"> -0.4608541 </td>
   <td style="text-align:right;"> -1.2596689 </td>
   <td style="text-align:right;"> -0.8004630 </td>
   <td style="text-align:right;"> -0.8584401 </td>
   <td style="text-align:right;"> -0.4918316 </td>
   <td style="text-align:right;"> -0.2745069 </td>
   <td style="text-align:right;"> 0.2867235 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NP_001001514.1@\hphantom{}K.ET\*RTSSES\*IVSVPASSTSGSPSR.V </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -0.5242609 </td>
   <td style="text-align:right;"> -0.2564000 </td>
   <td style="text-align:right;"> -0.1703421 </td>
   <td style="text-align:right;"> -0.8016035 </td>
   <td style="text-align:right;"> -0.6177977 </td>
   <td style="text-align:right;"> -1.8986707 </td>
   <td style="text-align:right;"> -0.9217045 </td>
   <td style="text-align:right;"> -1.1410773 </td>
   <td style="text-align:right;"> -0.4219050 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NP_001001514.1@\hphantom{}K.ETRTSS\*ESIVSVPASSTSGSPSR.V </td>
   <td style="text-align:right;"> -0.3952598 </td>
   <td style="text-align:right;"> -0.3773814 </td>
   <td style="text-align:right;"> -0.1113785 </td>
   <td style="text-align:right;"> -0.2246094 </td>
   <td style="text-align:right;"> -0.5676749 </td>
   <td style="text-align:right;"> -0.7561756 </td>
   <td style="text-align:right;"> -0.9173361 </td>
   <td style="text-align:right;"> -1.0505509 </td>
   <td style="text-align:right;"> -1.4130691 </td>
   <td style="text-align:right;"> -0.2275133 </td>
   <td style="text-align:right;"> -0.8394326 </td>
   <td style="text-align:right;"> -0.3245755 </td>
   <td style="text-align:right;"> -0.9394990 </td>
   <td style="text-align:right;"> -0.5237369 </td>
   <td style="text-align:right;"> -0.8749379 </td>
   <td style="text-align:right;"> -0.7834843 </td>
   <td style="text-align:right;"> -0.5268181 </td>
   <td style="text-align:right;"> -0.2208240 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NP_001001514.1@\hphantom{}K.ETRTSSES\*IVS\*VPASSTSGSPSR.V </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -0.1167323 </td>
   <td style="text-align:right;"> 0.0249531 </td>
   <td style="text-align:right;"> 0.1045728 </td>
   <td style="text-align:right;"> -0.6020616 </td>
   <td style="text-align:right;"> -0.1762188 </td>
   <td style="text-align:right;"> -0.3440353 </td>
   <td style="text-align:right;"> -0.2297419 </td>
   <td style="text-align:right;"> 0.1989361 </td>
   <td style="text-align:right;"> 0.4882626 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NP_001001514.1@\hphantom{}K.GDADT\*RTNSPDLDSQS\*LSLSSGADQEPLQR.M </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -0.0742108 </td>
   <td style="text-align:right;"> -0.2566411 </td>
   <td style="text-align:right;"> -0.4073881 </td>
   <td style="text-align:right;"> -1.1553706 </td>
   <td style="text-align:right;"> -0.7643052 </td>
   <td style="text-align:right;"> -1.3815904 </td>
   <td style="text-align:right;"> -1.3611955 </td>
   <td style="text-align:right;"> -1.0650441 </td>
   <td style="text-align:right;"> -0.0341310 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NP_001001514.1@\hphantom{}K.GDADTRTNSPDLDS\*QS\*LSLSSGADQEPLQR.M </td>
   <td style="text-align:right;"> -0.8822947 </td>
   <td style="text-align:right;"> -0.7983921 </td>
   <td style="text-align:right;"> -0.6040902 </td>
   <td style="text-align:right;"> -1.0568789 </td>
   <td style="text-align:right;"> -0.5751194 </td>
   <td style="text-align:right;"> -0.4792532 </td>
   <td style="text-align:right;"> -1.1417354 </td>
   <td style="text-align:right;"> -0.9907350 </td>
   <td style="text-align:right;"> -0.8330389 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
  </tr>
</tbody>
</table></div>

</br>

We will save the cross-tab for later sections.


```r
# Modify cross-tab for saving
crosstab <- crosstab %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("phospho_peptide")

# Save cross-tab
write.table(crosstab, file = "data/phospho_quant_crosstab.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
```

