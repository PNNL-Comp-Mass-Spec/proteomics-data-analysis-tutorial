# Isobaric Quantification: Phosphoproteomics {#iso-phospho}

This pipeline shows how to process phosphoproteomics TMT data with `PlexedPiper`, though it can be used for any type of post-translational modification (PTM) TMT data. We will use data package 3626, which is "PlexedPiperTestData phospho". In addition to PlexedPiper, we will also need MSnID (the basis for PlexedPiper), PNNL.DMS.utils to interface with PNNL's DMS, and Biostrings to create an `AAStringSet` object from a FASTA file. Since a lot of these steps are the same as in Section \@ref(iso-global), a lot of the details will be omitted.




```r
## Uncomment to install missing packages
# install.packages("remotes")
# library(remotes)
# install_github("PNNL-Comp-Mass-Spec/MSnID@pnnl-master")
# install_github("PNNL-Comp-Mass-Spec/PlexedPiper")
# install_github("PNNL-Comp-Mass-Spec/PNNL.DMS.utils")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("Biostrings")
library(MSnID)
library(PlexedPiper)
library(PNNL.DMS.utils)
library(Biostrings)
```




## Prepare MS/MS Identifications

### Read MS-GF+ Data


```r
# Read MS-GF+ data
data_package_num <- 3626 # phospho
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


### Correct Isotope Selection Error 


```r
# Correct for isotope selection error
msnid <- correct_peak_selection(msnid)
```


### Remove Unmodified Peptides

Generally, we will remove unmodified peptides before any sort of filtering steps; however, unmodified peptides will be removed automatically in Section \@ref(map-mod-sites), so this step can be skipped if we need to tally the number of modified and unmodified peptides toward the end of processing.

In this case, the phosphorylation of an amino acid is marked by a `*` appearing next in the sequence. We can filter out peptides that do not contain this symbol with `apply_filter`. In regular expressions, the `*` is a special character called a metacharacter that must be escaped with backslashes, and the backslashes must also be escaped, since they are enclosed within a nested string (`"''"`). For non-metacharacters, it is not necessary to include the backslashes.


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


### Remove Contaminants 


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


### Improve Phosphosite Localization

Phospho datasets involve AScore jobs for improving phosphosite localization. There should be one AScore job per data package. If the AScore job does not exist, see <a href="https://prismwiki.pnl.gov/wiki/AScore_Job_Creation">AScore Job Creation</a> for how to set it up. The fetched object is a data.frame that links datasets, scans and original PTM localization to newly suggested locations. Importantly, it contains `AScore` column that signifies the confidence of PTM assignment. AScore > 17 is considered confident.


```r
# Filter PTMs by Ascore - only for phospho data
ascore <- get_AScore_results(data_package_num)
msnid <- best_PTM_location_by_ascore(msnid, ascore)
```

### MS/MS ID Filter: Peptide Level


```r
# 1% FDR filter at the peptide level
msnid <- filter_msgf_data(msnid, level = "peptide", fdr.max = 0.01)
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

### MS/MS ID Filter: Protein Level


```r
# Get path to FASTA file
path_to_FASTA <- path_to_FASTA_used_by_DMS(data_package_num)

# Compute number of peptides per 1000 amino acids
msnid <- compute_num_peptides_per_1000aa(msnid, path_to_FASTA)

# 1% FDR filter at the protein level
msnid <- filter_msgf_data(msnid, level = "accession", fdr.max = 0.01)
show(msnid)
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  23 
## #PSMs: 72481 at 0.12 % FDR
## #peptides: 21266 at 0.26 % FDR
## #accessions: 9424 at 0.98 % FDR
```

### Inference of Parsimonious Protein Set

If a protein was detected in the global proteomics results, we may be more confident that it will appear in the PTM results. We can perform prioritized inference of the protein set to ensure that, if a protein is reported in the global cross-tab, and it is present in the PTM MSnID after filtering, it will be included in the final PTM MSnID. We set the proteins from the global cross-tab as the prior.


```r
# Proteins from global proteomics cross-tab
load("data/3442_global_proteins.RData")

# Prioritized inference of parsimonious protein set
msnid <- infer_parsimonious_accessions(msnid, unique_only = FALSE,
                                       prior = global_proteins)
show(msnid)
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  23 
## #PSMs: 72481 at 0.12 % FDR
## #peptides: 21266 at 0.26 % FDR
## #accessions: 2895 at 1.6 % FDR
```

### Map Sites to Protein Sequences {#map-mod-sites}

`MSnID::map_mod_sites` creates a number of columns describing mapping of the modification sites onto the protein sequences. The most important for the user is `SiteID`. `names(fst)` must match `accessions(msnid)`; usually, we will have to modify names to select everything before the first space.


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
names(fst) <- gsub(" .*$", "", names(fst))
# First 6 names
head(names(fst))
```

```
## [1] "NP_783171.2"    "NP_001101862.2" "NP_113721.4"    "NP_714948.1"   
## [5] "NP_001000704.1" "NP_001000638.1"
```

The names are in the proper format, so we can continue with the main mapping call.


```r
# Main mapping call
msnid <- map_mod_sites(object = msnid, 
                       fasta = fst, 
                       accession_col = "accession", 
                       peptide_mod_col = "peptide",
                       mod_char = "*", # modification character
                       site_delimiter = ";") # ; between multiple sites
```

Table \@ref(tab:phospho-msnid-table) shows the first 6 rows of the processed MS-GF+ output.

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:90%; "><table class="table table-hover table-condensed" style="font-size: 12px; width: auto !important; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">(\#tab:phospho-msnid-table)<left>First 6 rows of the processed MS-GF+ results.</left>
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
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
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
   <td style="text-align:left;"> NP_001071138.1-S14;T19 </td>
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
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
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
   <td style="text-align:left;"> XP_006232986.1-T250 </td>
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
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
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
   <td style="text-align:left;"> XP_006232986.1-T250 </td>
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
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
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
   <td style="text-align:left;"> XP_006232986.1-T250 </td>
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
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
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
   <td style="text-align:left;"> NP_112621.1-T253 </td>
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
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
   <td style="text-align:right;"> 0 </td>
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
   <td style="text-align:left;"> NP_112621.1-T253 </td>
  </tr>
</tbody>
</table></div>

</br>

### Remove Decoy PSMs


```r
# Remove Decoy PSMs
msnid <- apply_filter(msnid, "!isDecoy")
show(msnid)
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  23 
## #PSMs: 72391 at 0 % FDR
## #peptides: 21211 at 0 % FDR
## #accessions: 2849 at 0 % FDR
```

## Prepare Reporter Ion Intensities

### Read MASIC Output 


```r
# Read MASIC data
masic_data <- read_masic_data_from_DMS(data_package_num, 
                                       interference_score = TRUE)
```

### Filter MASIC Data 


```r
# Filter MASIC data
masic_data <- filter_masic_data(masic_data, 
                                interference_score_threshold = 0.5,
                                s2n_threshold = 0)
```

## Create Study Design Tables

Aside from the fractions table, the other study design tables can be the same as those created for the global proteomics data.


```r
# Create fractions table
datasets <- unique(intersect(msnid$Dataset, masic_data$Dataset))
fractions <- data.frame(Dataset = datasets) %>% 
  mutate(PlexID = gsub(".*_P_(S\\d{1})_.*", "\\1", Dataset))

# Use global samples and references tables
samples <- read.delim("data/3442_samples.txt")
references <- read.delim("data/3442_references.txt")
```


```r
# Save phospho fractions table
write.table(fractions, file = "data/3626_fractions.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
```


## Create Quantitative Cross-tab


```r
# Create cross-tab - aggregate to SiteID level
crosstab <- create_crosstab(msnid = msnid, 
                            reporter_intensities = masic_data,
                            aggregation_level = "SiteID",
                            fractions = fractions, 
                            samples = samples, 
                            references = references)
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:90%; "><table class="table table-hover table-condensed" style="font-size: 12px; margin-left: auto; margin-right: auto;">
<caption style="font-size: initial !important;">(\#tab:unnamed-chunk-21)<left>First 6 rows of the phospho quantitative cross-tab.</left>
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
   <td style="text-align:left;"> NP_001000283.1-Y132;S137 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0.0038421 </td>
   <td style="text-align:right;"> -0.7416622 </td>
   <td style="text-align:right;"> -0.2425074 </td>
   <td style="text-align:right;"> -0.4879016 </td>
   <td style="text-align:right;"> -0.7497171 </td>
   <td style="text-align:right;"> -0.5428785 </td>
   <td style="text-align:right;"> -0.7653521 </td>
   <td style="text-align:right;"> -0.1239261 </td>
   <td style="text-align:right;"> -0.2119856 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NP_001001064.1-Y129 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0.1594274 </td>
   <td style="text-align:right;"> 0.1932351 </td>
   <td style="text-align:right;"> -0.6587720 </td>
   <td style="text-align:right;"> -0.4461276 </td>
   <td style="text-align:right;"> -0.6370354 </td>
   <td style="text-align:right;"> -0.2976259 </td>
   <td style="text-align:right;"> -1.0408865 </td>
   <td style="text-align:right;"> -0.1385762 </td>
   <td style="text-align:right;"> -0.7301538 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NP_001001512.2-S241 </td>
   <td style="text-align:right;"> -0.5441083 </td>
   <td style="text-align:right;"> -0.5540885 </td>
   <td style="text-align:right;"> -0.0291611 </td>
   <td style="text-align:right;"> -0.5675833 </td>
   <td style="text-align:right;"> -0.1812456 </td>
   <td style="text-align:right;"> -0.803056 </td>
   <td style="text-align:right;"> -0.0427948 </td>
   <td style="text-align:right;"> -1.0398773 </td>
   <td style="text-align:right;"> -0.5236890 </td>
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
   <td style="text-align:left;"> NP_001001512.2-S242 </td>
   <td style="text-align:right;"> -0.2806018 </td>
   <td style="text-align:right;"> -0.3923442 </td>
   <td style="text-align:right;"> -0.2939827 </td>
   <td style="text-align:right;"> -0.6512914 </td>
   <td style="text-align:right;"> -0.3242804 </td>
   <td style="text-align:right;"> -1.210593 </td>
   <td style="text-align:right;"> -0.4914234 </td>
   <td style="text-align:right;"> -0.8869142 </td>
   <td style="text-align:right;"> -0.6421375 </td>
   <td style="text-align:right;"> -0.3079991 </td>
   <td style="text-align:right;"> -0.4801256 </td>
   <td style="text-align:right;"> -0.8072458 </td>
   <td style="text-align:right;"> -0.7588038 </td>
   <td style="text-align:right;"> -0.4458455 </td>
   <td style="text-align:right;"> -0.6640721 </td>
   <td style="text-align:right;"> -1.3958746 </td>
   <td style="text-align:right;"> -0.3091972 </td>
   <td style="text-align:right;"> -0.7046127 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NP_001001512.2-S699 </td>
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
   <td style="text-align:right;"> -0.4608541 </td>
   <td style="text-align:right;"> -1.2173509 </td>
   <td style="text-align:right;"> -0.8004630 </td>
   <td style="text-align:right;"> -1.2596689 </td>
   <td style="text-align:right;"> -0.4918316 </td>
   <td style="text-align:right;"> -0.8584401 </td>
   <td style="text-align:right;"> 0.2867235 </td>
   <td style="text-align:right;"> -0.2745069 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> NP_001001512.2-S746 </td>
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
   <td style="text-align:right;"> -1.6828088 </td>
   <td style="text-align:right;"> -2.2975481 </td>
   <td style="text-align:right;"> -1.8632844 </td>
   <td style="text-align:right;"> -2.3039510 </td>
   <td style="text-align:right;"> -1.5981389 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> -0.4736739 </td>
   <td style="text-align:right;"> -2.4732843 </td>
  </tr>
</tbody>
</table></div>

</br>


```r
# Save cross-tab
write.table(crosstab, file = "data/3662_phospho_crosstab.txt",
            sep = "\t", quote = FALSE, row.names = TRUE)
```


