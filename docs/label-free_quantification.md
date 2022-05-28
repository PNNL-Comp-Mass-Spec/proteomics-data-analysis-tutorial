# Spectral Counting

This is a generic spectral counting script for MS-GF+ Human/UniProt searches. Only modify the lines that change the data package number and the name of the final .xlsx file that will be saved, unless you know what you are doing.


```r
## Uncomment to install missing packages
# install.packages("devtools")
# library(devtools)
# install_github("PNNL-Comp-Mass-Spec/MSnID@pnnl-master")
# install_github("PNNL-Comp-Mass-Spec/PlexedPiper")
# install_github("PNNL-Comp-Mass-Spec/PNNL.DMS.utils")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("MSnbase")
# install.packages("writexl")
# install.packages("dplyr")
# install.packages("tibble")
library(MSnID)
library(PlexedPiper)
library(PNNL.DMS.utils)
library(MSnbase)
library(writexl)
library(dplyr)
library(tibble)
```



```r
# Data package number
data_package_num <- 3987
# Name of the final file to save
file_name <- "data/3987_spectral_counts.xlsx"
```

Do not modify anything below unless you know what you are doing.


```r
# Read MS-GF+ results from the DMS
m <- read_msgf_data_from_DMS(data_package_num = data_package_num)

# Filter to 1% FDR at the peptide level
m <- filter_msgf_data(m, level = "peptide", fdr.max = 0.01)

# UniProt to gene symbol conversion table
conv_tab <- fetch_conversion_table(organism_name = "Homo sapiens", 
                                   from = "UNIPROT", to = "SYMBOL")
```

When running `fetch_conversion_table`, if a prompt appears that requires an answer, type `yes` and press enter.


```r
# Modify accessions column of psms to use gene symbols
m <- remap_accessions(m, conv_tab, "\\|([^|-]+)(-\\d+)?\\|")

# Do the same remapping to the FASTA file
fst_path <- path_to_FASTA_used_by_DMS(data_package_num = data_package_num)
fst_path_2 <- remap_fasta_entry_names(
  path_to_FASTA = fst_path, conversion_table = conv_tab, 
  extraction_pttrn = "\\|([^|-]+)(-\\d+)?\\|"
)

# Compute the number of amino acids per 1000 and use that to filter
# to 1% FDR at the protein level
m <- compute_num_peptides_per_1000aa(m, fst_path_2)
m <- filter_msgf_data(m, "accession", fdr.max = 0.01)

# Parsimonious protein inference
m <- infer_parsimonious_accessions(m)
show(m) # Assessment of filtering quality
```

```
## MSnID object
## Working directory: "."
## #Spectrum Files:  8 
## #PSMs: 114125 at 0.012 % FDR
## #peptides: 23209 at 0.043 % FDR
## #accessions: 2362 at 0.38 % FDR
```

The results look reasonable, so we will continue on to spectral counting.


```r
# Remove decoys
m <- apply_filter(m, "!isDecoy")

# Convert m to an MSnSet
msnset <- as(m, "MSnSet")

# Spectral counting:
# Within each accession group, sum the values within columns.
msnset <- combineFeatures(msnset,
                          fData(msnset)$accession,
                          redundancy.handler = "multiple",
                          method = "sum",
                          cv = FALSE)

# Sort features from most to least abundant
tot_count <- rowSums(exprs(msnset))
msnset <- msnset[order(-tot_count), ]
```


```r
# Save exprs as an .xlsx file
msnset %>%
  exprs() %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  write_xlsx(path = file_name)
```

