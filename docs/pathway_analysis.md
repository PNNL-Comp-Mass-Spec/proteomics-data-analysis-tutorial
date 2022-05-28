# Pathway Analysis {#pathway-analysis}



In Section \@ref(DEA), we covered analysis at the individual feature level (protein, peptide, phosphoprotein, etc.). While this is useful, it is not without its own set of shortcomings. For instance, there may be no features that pass the significance threshold after correcting for multiple hypothesis testing. Alternatively, there may be many features that are statistically significant, and interpreting this list can be tedious and "prone to investigator bias toward a hypothesis of interest" [@maleki_gene_2020]. Another issue is that differential analysis fails to detect subtle, yet coordinated changes in groups of related features [@subramanian_gene_2005].

In order to address these, and other, issues, pathway analysis instead examines *a priori* defined **gene sets**---groups of genes that participate in the same biological pathway, share the same cellular location, etc. In this section, we will explore some common annotation databases, as well as two pathway analysis methods: Over-Representation Analysis (ORA) and Gene Set Enrichment Analysis (GSEA).

## Annotation Databases {#annotation-databases}

In this section, we will explore a few of the common annotation databases used for pathway analysis.

### Gene Ontology {#gene-ontology}

The Gene Ontology (GO) database is divided into three separate domains: Biological Process, Cellular Component, and Molecular Function (see the <a href="http://geneontology.org/docs/ontology-documentation/" title = "Gene Ontology overview">Gene Ontology overview</a> for more details regarding each domain). Each domain is structured as a directed acyclic graph (DAG) where nodes are terms and edges are the <a href="http://geneontology.org/docs/ontology-relations/#:~:text=Main%20relations%20used%20in%20GO" title = "Main relations used in GO">relations</a> between the terms (part of, is a, has part, regulates). Nodes can be connected to multiple child and parent nodes, where the group of genes annotated to a child node is a subset of those that are annotated to its parent node(s) [@noauthor_relations_2021; @goeman_multiple_2008].

#### Semantic Similarity {#semantic-similarity}

Due to the DAG structure of each domain, there is often redundancy in pathway analysis results. For example, suppose terms <a href="https://www.ebi.ac.uk/QuickGO/term/GO:0006119" title = "oxidative phosphorylation">GO:0006119</a>, <a href="https://www.ebi.ac.uk/QuickGO/term/GO:0009060" title = "aerobic respiration">GO:0009060</a>, and <a href="https://www.ebi.ac.uk/QuickGO/term/GO:0046034" title = "ATP metabolic process">GO:0046034</a> are significantly over-represented biological processes. <a href="https://www.ebi.ac.uk/QuickGO/term/GO:0009060" title = "aerobic respiration">GO:0009060</a> and <a href="https://www.ebi.ac.uk/QuickGO/term/GO:0046034" title = "ATP metabolic process">GO:0046034</a> are the parent terms of <a href="https://www.ebi.ac.uk/QuickGO/term/GO:0006119" title = "oxidative phosphorylation">GO:0006119</a>. Due to this relationship, the terms likely provide much of the same information, so the inclusion of all three terms in the output is unnecessary. In order to resolve this redundancy, we can calculate the **semantic similarity** between pairs of GO terms, which "assesses the likeness in meaning of two concepts" [@pesquita_semantic_2017]. Basically, if two terms are highly related, we can use some other criteria (such as adjusted p-value or level in the DAG) to retain only one of the terms. Below, we use the `GOSemSim` package to calculate the semantic similarity between the terms.


```r
## Calculate semantic similarity between GO terms
library(GOSemSim)
library(org.Hs.eg.db)

# GO DATA for measuring semantic similarity.
# keytype is "ENTREZID" by default and 
# information content is calculated (computeIC = TRUE)
semData <- godata(OrgDb = "org.Hs.eg.db", ont = "BP")
terms <- c("GO:0006119", "GO:0009060", "GO:0046034")
# measure = "Rel" is the default for clusterProfiler::simplify
# See code for clusterProfiler:::simplify_internal
sim <- mgoSim(GO1 = terms, GO2 = terms, semData = semData, 
              measure = "Rel", combine = NULL) 
```



If `measure` is `"Lin"`, `"Jiang"`, or `"Wang"`, the semantic similarity of a term with itself will be 1. This is not true for the other methods.

We can see from Table \@ref(tab:sem-sim-table) that <a href="https://www.ebi.ac.uk/QuickGO/term/GO:0009060" title = "aerobic respiration">GO:0009060</a> and <a href="https://www.ebi.ac.uk/QuickGO/term/GO:0046034" title = "ATP metabolic process">GO:0046034</a> have low semantic similarity, while <a href="https://www.ebi.ac.uk/QuickGO/term/GO:0006119" title = "oxidative phosphorylation">GO:0006119</a> is highly similar to its parent terms. This makes sense because the parent terms are not related/connected in the DAG.

Now that we have the semantic similarities, we can remove redundant terms. `clusterProfiler` has a function called `simplify` that will calculate semantic similarity and remove terms. By default, if there are two terms with a Wang semantic similarity greater than 0.7, `simplify` retains the term with the lowest adjusted p-value. See <a href="https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/" title = "use simplify to remove redundancy of enriched GO terms">this post</a> by Guangchuang Yu for more details on `clusterProfiler::simplify`.

#### GO Subsets/Slims {#go-slim}

Another way to handle the redundancy of GO terms is to use a <a href="http://geneontology.org/docs/go-subset-guide/" title = "Guide to GO subsets">GO slim</a>, which is a subset of more general or research-relevant terms from the GO. GO slims can be <a href="http://geneontology.org/docs/download-ontology/#subsets" title = "Download the ontology">downloaded</a> or the `biomaRt` package can be used to access GO slim accessions.


```r
## Create human GO slim
library(biomaRt)
library(clusterProfiler) # gcSample data
library(dplyr)

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                dataset = "hsapiens_gene_ensembl")
# Uncomment to determine which attributes to select in getBM()
# View(listAttributes(mart)) 

# The GO slim columns are goslim_goa_accession and goslim_goa_description.
# We will map from the Entrez IDs in gcSample to these attributes.
data(gcSample)
universe <- unique(unlist(gcSample))
GO_slim <- getBM(filters = "entrezgene_id",
                 attributes = c("entrezgene_id",
                                "goslim_goa_accession",
                                "goslim_goa_description"),
                 values = universe, # Subset to these Entrez IDs
                 mart = mart) %>% 
  # Convert entrezgene_id from integer to character
  mutate_all(as.character)
```







Unfortunately, not every GO accession maps to a domain when we use `biomaRt` (unsure why this is the case), so we won't be able to separate the terms. However, there are two ways that we can still use these GO slim accessions. Either follow the steps for using `fgsea::fora` (Section \@ref(ora-examples)) with a gene set list that has been subset to the GO slim accessions, or remove any non GO slim accessions from the final results and readjust the remaining p-values.

### Reactome

[Home - Reactome Pathway Database](https://reactome.org/)

Reactome is a manually-curated database for biological processes and pathways. As of version 80 (April 2021), it contains data on 15 species (notably, H. sapiens, M. musculus, and R. norvegicus). H. sapiens is the most highly-annotated organism with 2580 pathways.

### KEGG

The Kyoto Encyclopedia of Genes and Genomes (KEGG)

[KEGG: Kyoto Encyclopedia of Genes and Genomes](https://www.genome.jp/kegg/)

```{=html}
<!---
### Pfam

[Pfam: Home page](http://pfam.xfam.org/)
[Pfam Documentation](https://pfam-docs.readthedocs.io/en/latest/)
[profile Hidden Markov Models](https://watermark.silverchair.com/140755.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAvYwggLyBgkqhkiG9w0BBwagggLjMIIC3wIBADCCAtgGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQM8AYKH-bv0m2nHfTQAgEQgIICqYjz7XAVxAREVHGvJofH3peUmhTifphKa4BaiK7Xw4XUV41B2Eolc7IBL6Q9Oh0MAVCD0yHlxqBSwJ3cETaSRp60T6Bk9Cl6nUyjTUpq-4f8YW1pUgOxvd7OCLkE11Twa60Y9T6JFyIUuCxHs6UgbRBAbwNZ-6HG9-WRHxiWzCSf6SbrH74ouzlsabLh1OZMSoJ3dfbNrBYRz-bkonlmHUBASLKee1YA6Eg4EavR1Qa07g5FwDJhGVUgDw4IyhkyGAucfJTwl5_FxFGrs8-FEQ01LZqQDzZLsL0zcH7NsD8dB6PaWcP-lrH0TQpRK-WHhSC2WmL6rEoV_zwzUE9faUkLRzbg82siXNFUhgR9cwpe_ycHF-ffmIVMjpPWYvmMzrk4YfMfp_A0OsNbu0S7vSL7gqJ-uD2HyeJ8XkAZ1HXDjBKGHyQmQX6otAHfdNHbHUMl8GaFmgUPDkfj9swuUeQ0udhKRRssJNF8cqrs-xDzNQnB3-cGhi0bfV7PT1CZxGBMhlC4UsTu0hLZUfMxozxXfDjLY0kZEu66lTJL7CYGru4JCx7qFHvaQZg3icjHLOAiwov7v9-CimdXepQrZxZQ588N2ZlgsZ5f0xiZ5qdGkGwDyhdFmLPj2f_RFXRe5TuCWLU1DTkX5NyIUYGwGJECn6n_qSsA37n73ek--_eBPV-F6GBYQe1oRWP2S5SQTxP0r3yPyh7dIlqcqQxWozlsew1fGCMjRoUv1e3SFU7LZ0j8bUSXXEAgOxjjOBSxWu1QzgjLa-MRV0IQ61tDXfGdJhSuuZSxRwlxTbpZKQWfaS3s1jFU9Hm-B-9CvhmJQgb7elOeD0cbpVMeM03ebaK-NgPA0dX5slcIe1--BtOdKoKjepE7--RS4Z4Bls4ZI0d0valj-zsHiw)
--->
```

### MSigDB {#msigdb}

The Molecular Signatures Database (MSigDB) is a comprehensive resource of manually-curated gene sets divided into nine collections, as of [v7.5.1](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/MSigDB_v7.5.1_Release_Notes "MSigDB v7.5.1 Release Notes") [@liberzon_molecular_2015]. **Importantly, it contains non-redundant versions of the most up-to-date Gene Ontology, Reactome, and KEGG databases.** The method for eliminating term redundancy is described in the [v7.0 Release Notes](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/MSigDB_v7.0_Release_Notes#C5_.28Gene_Ontology_collection.29_-_Major_overhaul):

> We computed Jaccard coefficients for each pair of sets, and marked a pair as highly similar if its Jaccard coefficient was greater than 0.85. We then clustered highly similar sets into "chunks" using the `hclust` function from the R `stats` package according to their GO terms and applied two rounds of filtering for every "chunk". First, we kept the largest set in the "chunk" and discarded the smaller sets. This left "chunks" of highly similar sets of identical sizes, which we further pruned by preferentially keeping the more general set (i.e., the set closest to the root of the GO ontology tree).

**Note:** the Jaccard coefficient is used to measure the similarity of two sets, $A$ and $B$, and is calculated as $J(A, B) = \frac{|A \cap B|}{|A \cup B|}$.

## Over-Representation Analysis {#ora}

### Overview {#ora-overview}

[Over-representation analysis (ORA) is used to determine which *a priori* defined gene sets are more present (over-represented) in a subset of "interesting" genes than what would be expected by chance]{style="color:#3366ff"} [@huang_bioinformatics_2009].

For example, if 10% of all genes being considered are "interesting" (statistically different between conditions, clustered together, etc.), we expect that about 10% of every gene set will be in this "interesting" group of genes. If any gene sets are significantly more present in this group, we would say they are over-represented. Note that **we are not limited to analyzing genes**, though that is how most databases are set up. So long as we have a feature-to-feature-set map, we can perform ORA.

**I recommend ORA only when GSEA is not appropriate**, such as for analyzing clusters (from k-means clustering, WGCNA, Mclust, etc.), where the only information available are the group designations. See Section \@ref(ora-drawbacks) (ORA Drawbacks) for details.

### Mathematical Details {#ora-math}

For each gene set, an enrichment p-value is calculated using the Binomial distribution, Hypergeometric distribution, the Fisher exact test, or the Chi-square test. Although this list is not all-encompassing, these are the most popular statistical methods [@huang_bioinformatics_2009]. Below is the formula for calculating the enrichment p-value for a particular gene set using the Hypergeometric distribution.

::: math
$$
P(X\geq x) = 1 - P(X \leq x-1) = 1 - \sum\limits_{i=0}^{x-1}\frac{\hphantom{}{M \choose i }{N - M \choose n-i}}{N \choose n}
$$
:::

In this equation, $N$ is the number of background genes, $n$ is the number of "interesting" genes, $M$ is the number of genes that are annotated to a particular gene set $S$, and $x$ is the number of "interesting" genes that are annotated to $S$. The numerator of the sum is the number of samples of $n$ genes that can be taken from a population of $N$ genes where exactly $i$ of the genes are annotated to $S$ and $n-i$ are not annotated to $S$. The denominator of the sum is the total number of samples of size $n$ that can be taken from a population of size $N$.

For example, suppose we have a list of 8000 genes, of which 400 are members of the same cluster $C$. Also suppose that 100 of the 8000 genes are annotated to a particular gene set $S$. Of these 100 genes, 20 are members of $C$. The probability that 20 or more (up to 100) genes annotated to $S$ are in cluster $C$ by chance is given by

::: math
$$
P(X\geq 20) = 1 - P(X \leq 19) = 1-\sum \limits_{i=0}^{19}\frac{\hphantom{}{100 \choose i}{8000 - 100 \choose 400-i}}{8000 \choose 400} = 7.88 \times 10^{-8}
$$
:::

That is, it is extremely unlikely that 20 of the 100 genes from this set are grouped in cluster $C$ by chance (at least, prior to adjustment for multiple comparisons). The code to calculate this p-value is


```r
phyper(q = 20 - 1, m = 400, n = 8000 - 400, k = 100, lower.tail = FALSE)
```

After a p-value has been calculated for each of the applicable gene sets, a multiple testing adjustment should be applied.

### Drawbacks {#ora-drawbacks}

ORA is not recommended as a follow-up to differential-expression analysis for the reasons below. Use GSEA instead, if appropriate.

1.  The choice of the threshold for statistical significance and the multiple comparison adjustment method can greatly impact the analysis [@huang_bioinformatics_2009].

2.  ORA fails to incorporate direction of change. (Are the genes in a given set mainly up or down-regulated in one condition relative to another?). It is NOT a good idea to split DEA results by the direction of change and apply ORA to the resulting subsets, unless you are specifically asking "which gene sets are over-represented when we only consider genes that are up- or down-regulated?"

3.  If few genes are in the "interesting" group, ORA may not yield useful or reliable results. For example, suppose 30 out of 8000 genes are "interesting". 100 of the genes are annotated to a particular gene set, of which 3 are "interesting". The associated Hypergeometric p-value is 0.006, and this set would be considered significantly over-represented at the 0.01 level (at least, prior to p-value adjustment); however, if only 2 of the genes in this set are "interesting", this p-value increases 10-fold to 0.0536 and is no longer significant even at the 0.05 level.

4.  ORA can not be used if the input contains duplicates. For example, a single feature can not be a member of two or more groups or present multiple times in the same group. This usually happens when attempting to perform gene-level ORA on protein-level differential analysis results, and can lead to artificial over-representation if genes are counted multiple times. Instead, use GSEA and summarize the ranking metric at the gene level (take the average).

### Examples {#ora-examples}

For these examples, we will show how to perform ORA with the `fgsea`, `clusterProfiler`, `ReactomePA`, and `GOstats` packages on clustering results. The databases that we will cover are Gene Ontology Biological Processes and Reactome, and we will only consider gene sets/pathways with at least 15 and no more than 300 genes. For details on these different annotation databases, please see Section \@ref(annotation-databases) (Annotation Databases).

We will use the `gcSample` data from `clusterProfiler` and treat the entire list as the **gene universe/background**. Each gene is represented by a human Entrez gene ID, which is the default keytype used by the `clusterProfiler` functions (and the only keytype compatible with `ReactomePA::enrichPathway`).


```r
library(clusterProfiler)
data("gcSample") # data for examples
```

Since the genes in `gcSample` are not unique, we will subset to unique genes for the sake of these examples.


```r
# Need to remove duplicates for the examples
all_genes <- unlist(gcSample)
universe <- all_genes[Biobase::isUnique(all_genes)] # all unique genes

# List with only unique genes
gcUnique <- lapply(gcSample, function(group_i) {
  group_i[group_i %in% universe]
})
```

#### ORA with `fgsea` {#fgsea-ora}

The `fgsea` package can be used to perform over-representation analysis with the `fora` function, which applies the hypergeometric test. It requires a list of gene sets or pathways, a vector of "interesting" genes to test, and a gene universe vector. For clustering results, we can run this function in `lapply` to test each cluster. For this example, we will perform ORA on the non-redundant Gene Ontology biological processes (GO-BP) sets from MSigDB, but any database can be used.

We first need to get the list of gene sets; we do this with the `msigdbr` package (Section \@ref(msigdb)). We use `msigdbr_collections` to determine the category and subcategory and `msigdbr` to fetch the data.


```r
# MSigDB R package
library(msigdbr)
msigdbr::msigdbr_collections() # available collections
# Subset to Human GO-BP sets
BP_db <- msigdbr(species = "Homo sapiens", 
                 category = "C5", subcategory = "GO:BP")
head(BP_db)
```


```r
# Convert to a list of gene sets
BP_conv <- unique(BP_db[, c("entrez_gene", "gs_exact_source")])
BP_list <- split(x = BP_conv$entrez_gene, f = BP_conv$gs_exact_source)
# First ~6 IDs of first 3 terms
lapply(head(BP_list, 3), head)
```

We have all the required input, so we can move on to ORA. In order to perform ORA on each cluster, we can wrap `fora` in `lapply`.


```r
## Cluster GO-BP ORA with fgsea package
library(fgsea)
library(dplyr)

# For each cluster i, perform ORA
fgsea_ora <- lapply(seq_along(gcUnique), function(i) {
  fora(pathways = BP_list, 
       genes = gcUnique[[i]], # genes in cluster i
       universe = universe, # all genes
       minSize = 15, 
       maxSize = 500) %>% 
    mutate(cluster = names(gcUnique)[i]) # add cluster column
}) %>% 
  data.table::rbindlist() %>% # combine tables
  filter(padj < 0.05) %>% 
  arrange(cluster, padj) %>% 
  # Add additional columns from BP_db
  left_join(distinct(BP_db, gs_subcat, gs_exact_source, 
                     gs_name, gs_description),
            by = c("pathway" = "gs_exact_source")) %>% 
  # Reformat descriptions
  mutate(gs_name = sub("^GOBP_", "", gs_name),
         gs_name = gsub("_", " ", gs_name))

# First 6 rows
head(fgsea_ora)
```

See `?fora` for a description of the output. The base output does not include term descriptions, so we had to add those ourselves.


#### ORA with `clusterProfiler`/`ReactomePA` {#ora-clustprof}

`clusterProfiler` and `ReactomePA` are convenience packages that use `fgsea` and `AnnotationDbi` as the basis for their pathway analysis functions. The functions are a bit more user-friendly, as the user does not need to fetch the gene lists themselves, but they are much slower, and the resulting objects are not simple data.frames.

`clusterProfiler` provides the `enrichGO` and `enrichKEGG` functions for GO and KEGG ORA, respectively (among others). `ReactomePA` provides the `enrichPathway` function for Reactome pathway ORA. For other databases, the `clusterProfiler::enricher` function can be used (though this is slower than just using `fgsea::fora` with a list of gene sets). To perform ORA on clustering results, we use `clusterProfiler::compareCluster` and tell it to use `enrichGO` as the ORA function.


```r
## Cluster GO-BP ORA with clusterProfiler package
cp_ora <- compareCluster(
  geneClusters = gcUnique, 
  fun = "enrichGO", # ORA function to apply to each cluster
  # Arguments below are passed to enrichGO
  OrgDb = "org.Hs.eg.db", 
  keyType = "ENTREZID", 
  ont = "BP", # BP, CC, MF, or ALL for all ontologies
  pvalueCutoff = 0.05,
  qvalueCutoff = 1, # do not filter by q-value
  pAdjustMethod = "BH", # p-values are adjusted within clusters
  universe = universe, # all genes
  minGSSize = 15, 
  maxGSSize = 500
)

# First 6 entries sorted by cluster and p-value
cp_ora@compareClusterResult %>% 
  arrange(Cluster, pvalue) %>% 
  head()
```

Unlike the `fgsea::fora` results, these include the description of each term. As for the other columns: `GeneRatio` is the same as `overlap` (from the `fora` results) divided by the cluster size, `BgRatio` is the set size divided by the universe size, `pvalue` is the raw p-value, `p.adjust` is the BH-adjusted p-value, `qvalue` is the q-value, `geneID` is the same as `overlapGenes` from `fora`, and `Count` is the overlap size.

Below is an example of how to perform Reactome ORA with `ReactomePA::enrichPathway`.


```r
## Reactome ORA with ReactomePA package
library(ReactomePA)
react_ora <- compareCluster(
  geneClusters = gcUnique, 
  fun = "enrichPathway", # ORA function to apply to each cluster
  # Arguments below are passed to enrichPathway
  organism = "human",
  pvalueCutoff = 1, # Do not filter by p-value
  qvalueCutoff = 1, # Do not filter by q-value
  pAdjustMethod = "BH", # p-values are adjusted within clusters
  universe = universe, # all genes
  minGSSize = 15, 
  maxGSSize = 500
)
# First 6 rows
head(react_ora@compareClusterResult)
```

#### ORA with `GOstats`

In the previous ORA examples, the Hypergeometric test is performed independently for each gene set; however, this does not capture the relationship between GO terms (described in Section \@ref(gene-ontology)). Since "each GO term inherits all annotations from its more specific descendants," results tend to be redundant (except when using MSigDB), as they include directly-related GO terms with a high degree of overlap [@falcon_using_2007]. One way to handle this is with a procedure that conditions on the GO structure, like the one described by S. Falcon and R. Gentleman [-@falcon_using_2007]:

> Given a subgraph of one of the three GO ontologies [BP, MF, or CC], we test the leaves of the graph, that is, those terms with no child terms. Before testing the terms whose children have already been tested, we remove all genes annotated at significant children [`pvalueCutoff = 0.05` in the code below] from the parent's gene list. This continues until all terms have been tested.

This approach is implemented in the `GOstats` package by setting `conditional = TRUE` when creating a new object of class `GOHyperGParams` (See `help("GOHyperGParams-class", package = "Category")` for more details). Below, we will perform GO-BP ORA just for the first cluster of `gcUnique` because it is time-consuming, but this could be wrapped in `lapply` to get results for each cluster (like in Section \@ref(fgsea-ora)).


```r
## Cluster GO-BP ORA with GOstats package
library(GOstats)
library(org.Hs.eg.db)
library(dplyr)
# For cluster 1, perform conditional ORA
gostats_ora <- new(Class = "GOHyperGParams",
                   ontology = "BP",
                   geneIds = gcUnique[[1]],
                   universeGeneIds = universe,
                   annotation = "org.Hs.eg.db",
                   pvalueCutoff = 0.05,
                   testDirection = "over",
                   conditional = TRUE, # condition on GO structure
                   minSizeCutoff = 15,
                   maxSizeCutoff = 500) %>% 
  hyperGTest() %>% # Hypergeometric testing
  summary() %>% # extract results
  # adjust p-values
  mutate(Padj = p.adjust(Pvalue, method = "BH")) %>% 
  filter(Padj < 0.05) %>% 
  arrange(Padj)
# First 6 rows
head(gostats_ora)
```

If we compare this table to the results from the other packages, we see that only <GO:0007600> and <GO:0061844> from the top 6 terms made it into the above table.

## Gene Set Enrichment Analysis {#gsea}

### Overview {#gsea-overview}

**(More details to be added at a later date.)**

[Gene set enrichment analysis (GSEA) is a rank-based approach that determines whether predefined groups of genes/proteins/etc. are primarily up or down in one condition relative to another]{style="color:#3366ff"} [@mootha_pgc-1-responsive_2003; @subramanian_gene_2005]. It is typically performed as a follow-up to differential analysis, and is preferred to ORA (Section \@ref(ora)).

### Examples {#gsea-examples}

These examples will show how to run Fast GSEA (FGSEA) in R, which is based on the gene permutation approach [@korotkevich_fast_2016].

The input of FGSEA is a list of gene sets/pathways to check and a uniquely-named vector of ranking metric values sorted in descending order. The ranking metric that we will use is $-log_{10}(\text{p-value}) * sign(\text{logFC})$, but we could have easily used t-statistics or some other metric. We will perform differential analysis on the `cptac_oca` data from the `MSnSet.utils` package and use the results to create the ranking metric vector. We start off by remapping features from RefSeq to Entrez ID.


```r
## Fetch REFSEQ to ENTREZID conversion table
library(MSnID)
conv_tbl <- fetch_conversion_table(organism_name = "Homo sapiens", 
                                   from = "REFSEQ", to = "ENTREZID")
head(conv_tbl)
```

Now, we will create the differential analysis results and add the ENTREZID column.


```r
library(MSnSet.utils)
data("cptac_oca")
m1 <- oca.set

# Differential analysis
res <- limma_a_b(eset = m1, 
                 model.str = "~ PLATINUM.STATUS",
                 coef.str = "PLATINUM.STATUS")
# table(res$adj.P.Val < 0.05) # 0
# hist(res$P.Value, breaks = seq(0, 1, 0.05)) # looks uniform
head(res)
```


```r
library(dplyr)
# Add ENTREZID column
res <- res %>% 
  mutate(REFSEQ = sub("\\.\\d+", "", rownames(.))) %>% 
  left_join(conv_tbl, by = "REFSEQ")
head(res)
```

We need one ranking metric value per Entrez ID, so we will calculate the ranking metric and then take the gene-wise average. Then, we sort the genes in descending order by their ranking metric values and convert to a named vector called `geneList`.


```r
## Ranking metric vector for GSEA
geneList <- res %>%
  filter(!is.na(ENTREZID), !is.na(logFC)) %>% 
  mutate(ranking_metric = -log10(P.Value)*sign(logFC)) %>% 
  group_by(ENTREZID) %>% 
  summarise(ranking_metric = mean(ranking_metric, na.rm = TRUE)) %>% 
  arrange(-ranking_metric) %>% # sort descending (important!)
  tibble::deframe() # convert to named vector
head(geneList)
tail(geneList)
```

Now that we have the uniquely-named ranking metric vector, we can proceed with FGSEA. For these examples, we will show how to perform FGSEA with the `fgsea`, `clusterProfiler`, and `ReactomePA` packages. The databases that we will cover are Gene Ontology Biological Processes and Reactome, and we will only consider gene sets/pathways with at least 15 and no more than 500 genes. The maximum gene set size, as well as the number of permutations, affect the normalized enrichment scores and p-values. For details on these different annotation databases, please see Section \@ref(annotation-databases) (Annotation Databases).

#### GSEA with `fgsea`

To perform FGSEA with the `fgsea` package, we need a list of gene sets/pathways and the ranking metric vector. Below is one way to get the gene set list.


```r
# MSigDB R package
library(msigdbr)
msigdbr::msigdbr_collections() # available collections
# Subset to Human GO-BP sets
BP_db <- msigdbr(species = "Homo sapiens", 
                 category = "C5", subcategory = "GO:BP")
head(BP_db)
```


```r
# Convert to a list of gene sets
BP_conv <- unique(BP_db[, c("entrez_gene", "gs_exact_source")])
BP_list <- split(x = BP_conv$entrez_gene, f = BP_conv$gs_exact_source)
# First ~6 IDs of first 3 terms
lapply(head(BP_list, 3), head)
```


The `fgseaMultilevel` function uses the adaptive multilevel split Monte Carlo approach described in the original FGSEA paper [@korotkevich_fast_2016].


```r
## GO-BP FGSEA with fgsea package
library(fgsea)
set.seed(99)
system.time( # keep track of elapsed time
  fgsea_res <- fgseaMultilevel(pathways = BP_list,
                               stats = geneList, 
                               minSize = 15, 
                               maxSize = 500, 
                               eps = 0, 
                               nPermSimple = 10000)
)
# First 6 rows with lowest enrichment p-values
fgsea_res %>% 
  # Add additional columns from BP_db
  left_join(distinct(BP_db, gs_subcat, gs_exact_source, 
                     gs_name, gs_description),
            by = c("pathway" = "gs_exact_source")) %>% 
  # Reformat descriptions
  mutate(gs_name = sub("^GOBP_", "", gs_name),
         gs_name = gsub("_", " ", gs_name)) %>% 
  arrange(padj) %>% 
  head(8)
```

See `?fgseaMultilevel` for a description of the output. The output does not include term descriptions, so we had to add those ourselves.

#### GSEA with `clusterProfiler`/`ReactomePA`

See Section \@ref(ora-clustprof) to get a better understanding of these packages. `clusterProfiler` provides the `gseGO` and `gseKEGG` functions (among others) for FGSEA of the GO and KEGG databases, respectively. They are essentially more user-friendly wrapper functions that make use of the `fgsea` and `AnnotationDbi` package, but they tend to be much slower.


```r
## GO-BP FGSEA with clusterProfiler package
library(clusterProfiler)
system.time( # keep track of elapsed time
  cgsea_res <- gseGO(geneList = geneList, 
                     ont = "BP", 
                     OrgDb = "org.Hs.eg.db", 
                     minGSSize = 15, 
                     maxGSSize = 500, 
                     eps = 0, 
                     nPermSimple = 10000, 
                     seed = TRUE)
)
# First 8 rows with lowest enrichment p-values
cgsea_res@result %>% 
  arrange(pvalue) %>% 
  head(8)
```

Notice that the normalized enrichment scores (NES) are not quite the same as what we got when we used `fgseaMultiLevel`. This has to do with the permutations being different. To increase the precision of the NES and p-values, we can increase `nPermSimple`, though 10,000 should be more than sufficient.

Now, here is GSEA with the Reactome database using the `gsePathway` function from the `ReactomePA` package.


```r
## Reactome FGSEA with ReactomePA package
library(ReactomePA)
fgsea_react <- gsePathway(geneList = geneList, 
                          organism = "human",
                          minGSSize = 15, 
                          maxGSSize = 500, 
                          eps = 0, 
                          nPermSimple = 10000, 
                          seed = TRUE)
# First 6 rows
head(fgsea_react@result)
```

