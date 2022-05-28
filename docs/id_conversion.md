# Feature ID Conversion {#id-conversion}

This section shows how to convert from one feature ID (UniProt accessions, gene symbols, etc.) to another using different packages.




## Conversion with `biomaRt`


```r
## Install missing packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
bio_pkgs <- c("biomaRt", "Biostrings")
for (pkg_i in bio_pkgs) {
  if (!require(pkg_i, quietly = T, character.only = T))
    BiocManager::install(pkg_i)
}
if (!require("remotes", quietly = T)) install.packages("remotes")
if (!require("MSnID", quietly = T, character.only = T))
  remotes::install_github("PNNL-Comp-Mass-Spec/MSnID@pnnl-master")
## ------------------------
library(biomaRt) # ID conversion
library(Biostrings) # read FASTA files
library(MSnID) # parse_FASTA_names
```

The first steps are to determine which mart and dataset to use. `listMarts` will show the available marts. The first 6 rows of the available datasets (provided by `listDatasets(mart)`) are also shown. (Use `View`, rather than `head`, to search for the desired database.)


```r
# Create mart
listMarts() # determine biomart for useMart
```

```
##                biomart                version
## 1 ENSEMBL_MART_ENSEMBL      Ensembl Genes 106
## 2   ENSEMBL_MART_MOUSE      Mouse strains 106
## 3     ENSEMBL_MART_SNP  Ensembl Variation 106
## 4 ENSEMBL_MART_FUNCGEN Ensembl Regulation 106
```

```r
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL")
head(listDatasets(mart)) # determine dataset for useMart
```

```
##                        dataset                           description
## 1 abrachyrhynchus_gene_ensembl Pink-footed goose genes (ASM259213v1)
## 2     acalliptera_gene_ensembl      Eastern happy genes (fAstCal1.2)
## 3   acarolinensis_gene_ensembl       Green anole genes (AnoCar2.0v2)
## 4    acchrysaetos_gene_ensembl       Golden eagle genes (bAquChr1.2)
## 5    acitrinellus_gene_ensembl        Midas cichlid genes (Midas_v5)
## 6    amelanoleuca_gene_ensembl       Giant panda genes (ASM200744v2)
##       version
## 1 ASM259213v1
## 2  fAstCal1.2
## 3 AnoCar2.0v2
## 4  bAquChr1.2
## 5    Midas_v5
## 6 ASM200744v2
```

```r
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                dataset = "rnorvegicus_gene_ensembl")
```

Next, we determine which attributes to select for the conversion table with `listAttributes(mart)`. From that table, we select "refseq_peptide" and "external_gene_name".


```r
# Create conversion table
head(listAttributes(mart)) # determine attributes for getBM
```

```
##                            name                  description         page
## 1               ensembl_gene_id               Gene stable ID feature_page
## 2       ensembl_gene_id_version       Gene stable ID version feature_page
## 3         ensembl_transcript_id         Transcript stable ID feature_page
## 4 ensembl_transcript_id_version Transcript stable ID version feature_page
## 5            ensembl_peptide_id            Protein stable ID feature_page
## 6    ensembl_peptide_id_version    Protein stable ID version feature_page
```

```r
conv_tbl1 <- getBM(attributes = c("refseq_peptide", "external_gene_name"),
                   mart = mart)
head(conv_tbl1, 10)
```

```
##    refseq_peptide external_gene_name
## 1                                   
## 2                         AC118165.1
## 3    NP_001000130              Olr56
## 4    NP_001000302             Olr473
## 5                         AC099294.1
## 6                     AABR07054368.1
## 7                             Olr760
## 8    NP_001014048              Clrn3
## 9                     AABR07000137.1
## 10   NP_001011937              Doc2g
```

This table has a lot of blank entries that need to be removed.


## Conversion with `AnnotationHub`


```r
## Install missing packages
if (!require("remotes", quietly = T)) install.packages("remotes")
if (!require("MSnID", quietly = T))
  remotes::install_github("PNNL-Comp-Mass-Spec/MSnID@pnnl-master")
```

`MSnID` has a function called `fetch_conversion_table` that leverages `AnnotationHub` to create a conversion table.


```r
# Create conversion table with MSnID::fetch_conversion_table
conv_tbl2 <- fetch_conversion_table(
  organism_name = "Rattus norvegicus", from = "REFSEQ", to = "SYMBOL"
)
head(conv_tbl2, 10)
```

```
##          REFSEQ SYMBOL
## 1     NP_443211   Asip
## 2     NP_036620    A2m
## 3  XP_038962922    A2m
## 4  XP_038962923    A2m
## 5  NP_001382591 Acaa1a
## 6     NP_036621 Acaa1a
## 7  XP_038936714 Acaa1a
## 8  XP_038936715 Acaa1a
## 9     NP_058682  Acadm
## 10 NP_001104565   Acly
```


## Conversion Using FASTA Headers

If specifically converting to gene symbols, it is recommended to use the information in the headers of the FASTA file that was used for the database search. The gene symbol is always given by `GN=...`, so we can use a regular expression to extract it. For UniProt FASTA files, there is a function in `MSnID` called `parse_FASTA_names` that will extract the components of the FASTA headers and create a `data.frame`.


```r
## Read FASTA file
fst_path <- system.file("extdata/uniprot_rat_small.fasta.gz",
                        package = "MSnID")
conv_tbl3 <- parse_FASTA_names(path_to_FASTA = fst_path)
head(conv_tbl3)
```

```
##               feature database uniprot_acc isoform entry_name
## 1  sp|P63088|PP1G_RAT       sp      P63088      NA   PP1G_RAT
## 2 sp|Q4FZV7|TMUB2_RAT       sp      Q4FZV7      NA  TMUB2_RAT
## 3 sp|O55159|EPCAM_RAT       sp      O55159      NA  EPCAM_RAT
## 4 sp|Q80VJ4|GPCP1_RAT       sp      Q80VJ4      NA  GPCP1_RAT
## 5 sp|Q66MI6|T10IP_RAT       sp      Q66MI6      NA  T10IP_RAT
## 6 sp|O70453|HMOX3_RAT       sp      O70453      NA  HMOX3_RAT
##                                                        description
## 1 Serine/threonine-protein phosphatase PP1-gamma catalytic subunit
## 2     Transmembrane and ubiquitin-like domain-containing protein 2
## 3                                Epithelial cell adhesion molecule
## 4                   Glycerophosphocholine phosphodiesterase GPCPD1
## 5                   Testis-specific protein 10-interacting protein
## 6                                        Putative heme oxygenase 3
##            organism organism_id     gene protein_existence sequence_version
## 1 Rattus norvegicus       10116   Ppp1cc                 1                1
## 2 Rattus norvegicus       10116    Tmub2                 2                1
## 3 Rattus norvegicus       10116    Epcam                 1                1
## 4 Rattus norvegicus       10116   Gpcpd1                 1                1
## 5 Rattus norvegicus       10116 Tsga10ip                 2                2
## 6 Rattus norvegicus       10116    Hmox3                 5                1
```

