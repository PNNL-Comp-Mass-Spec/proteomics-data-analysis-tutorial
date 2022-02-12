# Feature ID Conversion {#id-conversion}

This section shows how to convert from one feature ID (UniProt accessions, gene symbols, etc.) to another using different packages.




## Conversion with `biomaRt`


```r
## Install missing packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
## ------------------------
library(biomaRt)
```


```r
# Create mart
listMarts() # determine biomart for useMart
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL")
View(listDatasets(mart)) # determine dataset for useMart
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                dataset = "rnorvegicus_gene_ensembl")
```


```r
# Create conversion table
View(listAttributes(mart)) # determine attributes for getBM
conv_tbl1 <- getBM(attributes = c("refseq_peptide", "external_gene_name"),
                   mart = mart)
head(conv_tbl1, 10)
```

## Conversion with `AnnotationHub`


```r
## Install missing packages
if (!require("remotes", quietly = T)) install.packages("remotes")
if (!require("MSnID", quietly = T))
  remotes::install_github("PNNL-Comp-Mass-Spec/MSnID@pnnl-master")
```



```r
# Create conversion table with MSnID::fetch_conversion_table
conv_tbl2 <- fetch_conversion_table(
  organism_name = "Rattus norvegicus", from = "REFSEQ", to = "SYMBOL"
)
head(conv_tbl2, 10)
```

