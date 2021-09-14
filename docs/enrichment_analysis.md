# Enrichment Analysis {#enrich}

<!---
TODO
* Explain why I am talking about genes instead of proteins.
* Create a better hook for this chapter. The transition feels too abrupt.
--->

While single-gene analysis (Section \@ref(DEA)) is a useful tool, it is not without its own set of shortcomings. For instance, there may be no genes that pass the significance threshold after correcting for multiple hypothesis testing. Alternatively, there may be many genes that are statistically significant, and interpreting this list can be tedious and "prone to investigator bias toward a hypothesis of interest" [@maleki_gene_2020]. Another issue is that single-gene analysis fails "to detect biological processes...that are distributed across an entire network of genes and subtle at the level of individual genes." That is, it is liable to miss effects that are the result of small changes in many related genes. In order to address these and other issues, enrichment analysis (also called over-representation analysis) analyzes **gene sets**—"groups of genes that share [a] common biological function, chromosomal location," involvement in a pathway, etc.—rather than individual genes [@subramanian_gene_2005].

<!---
TODO:
* Talk about "General principle of enrichment analysis"
--->
@huang_bioinformatics_2009 classifies enrichment analysis algorithms into three groups: singular enrichment analysis (SEA), gene set enrichment analysis (GSEA), and modular enrichment analysis (MEA).


## SEA

Singular Enrichment Analysis (SEA) is used to determine which gene sets are over-represented in a set of "interesting" genes compared to a set of background genes. For each gene set, an enrichment p-value is calculated using the Binomial distribution, Hypergeometric distribution, the Fisher exact test, or the Chi-square test. Although this list is not all-encompassing, these are the most popular statistical methods [@huang_bioinformatics_2009]. Below is the formula for calculating the enrichment p-value for a particular gene set using the Hypergeometric distribution.

\[
P(X\geq x) = 1 - P(X \leq x-1) = 1 - \sum\limits_{i=0}^{x-1}\frac{\hphantom{}{M \choose i }{N - M \choose n-i}}{N \choose n}
\]

In this equation, $N$ is the number of background genes, $n$ is the number of "interesting" (i.e. statistically-significant) genes, $M$ is the number of genes that are annotated to a particular gene set $S$, and $x$ is the number of "interesting" genes that are annotated to $S$ (i.e. $x = M \bigcap n$). The numerator of the sum is the number of samples of $n$ genes that can be taken where exactly $i$ of the genes are annotated to $S$ and $n-i$ are not annotated to $S$. The denominator of the sum is the number of samples of size $n$ that can be taken from a population of size $N$.

For example, suppose we have a list of 8000 genes, of which 400 are differentially expressed. Also suppose that 100 of the 8000 genes are annotated to a particular gene set $S$. Of these 100 genes, 20 are differentially expressed. The probability that 20 or more (up to 100) genes annotated to $S$ are differentially expressed by chance is given by

\[
P(X\geq 20) = 1 - P(X \leq 19) = 1-\sum \limits_{i=0}^{19}\frac{\hphantom{}{100 \choose i}{8000 - 100 \choose 400-i}}{8000 \choose 400} = 7.88 \times 10^{-8}
\]

That is, it is unlikely that $S$ is enriched by chance. The code to calculate this p-value is


```r
phyper(q = 20 - 1, m = 400, n = 8000 - 400, k = 100, lower.tail = FALSE)
```

## GSEA

<!---
Gene Set Enrichment Analysis (GSEA) employs a "no-cutoff" strategy that utilizes some "quantitative biological value (such as fold change or degree of differential expression)" to rank a list $L$ of genes [@tipney_introduction_2010; @subramanian_gene_2005]. After ranking $L$, an enrichment score (ES) is calculated for a particular gene set $S$. This is done by 1) "walking down the list $L$, increasing a running-sum statistic when we encounter a gene in $S$ and decreasing it when we encounter genes not in $S$"; and then 2) determining the "maximum deviation from zero encountered in the random walk." This maximum deviation is the ES for a particular gene set. Once the ES has been obtained, the statistical significance is estimated by permuting (randomizing) the "phenotype labels and recomput[ing] the ES of the gene set for the permuted data, which generates a null distribution for the ES" [@subramanian_gene_2005]. Finally, adjustment for multiple hypothesis testing is performed.
--->


Gene Set Enrichment Analysis (GSEA) employs a "no-cutoff" strategy that utilizes some "quantitative biological value (such as fold change or degree of differential expression)" to rank a list of $N$ genes in descending order [@tipney_introduction_2010]. Using this ranked list $L$, the values of the ranking metric, and an *a priori* defined gene set $S$, we can calculate an enrichment score. This is done by "walking down the list $L$, increasing a running-sum statistic when we encounter a gene in $S$ and decreasing it when we encounter genes not in $S$." The "maximum deviation from zero" of this running-sum statistic is the enrichment score for $S$, denoted by $ES(S)$. If the genes in $S$ are randomly distributed throughout $L$ (this is the null hypothesis), $ES(S)$ will be relatively small; however, if they are not randomly distributed (i.e. primarily located near the top or bottom of $L$), then $ES(S)$ will be relatively large [@subramanian_gene_2005].

To assess the statistical significance of $ES(S)$, we randomize the phenotype labels in the original expression data, re-compute the ranking metric values, reorder $L$ using these new values, and re-compute the enrichment score $ES_{\text{NULL}}$. This process is repeated hundreds or thousands of times until we have a distribution of enrichment scores. We then define the nominal p-value as the proportion of these enrichment scores that are greater than or equal to $ES(S)$. Once a p-value has been calculated for each gene set, they are adjusted for multiply hypothesis testing [@subramanian_gene_2005].

<!---
We will now elaborate on how $ES(S)$ is derived. Given a ranked list of genes $L = \{g_1, \dots, g_N \}$, a gene set $S$ consisting of $M$ genes, and the values of the ranking metric ($X = \{ x_1, \dots, x_N \}$ used to generate $L$, we define

\[
P_{\text{hit}}(S, i) = \sum \limits_{g_j \in S \\ j \leq i} \frac{|x_j|^p}{N_R} \quad \text{and} \quad P_{\text{miss}}(S, i) = \sum \limits_{g_j \in S \\ j \leq i} \frac{1}{N - M}
\]

where $N_R = \sum \limits_{g_j \in S} |x_j|^p$
--->

<!---
TODO:
Talk about benefits of GSEA over SEA (i.e. reduced variability in "non-statistical layers" and usable when few or no genes pass the threshold for statistical significance.)
--->


## Databases

### Gene Ontology

The Gene Ontology (GO) database is divided into three separate ontologies: Biological Process, Cellular Component, and Molecular Function.

[Gene Ontology overview](http://geneontology.org/docs/ontology-documentation/)

### Reactome

[Home - Reactome Pathway Database](https://reactome.org/)

### KEGG

[KEGG: Kyoto Encyclopedia of Genes and Genomes](https://www.genome.jp/kegg/)

### Pfam

[Pfam: Home page](http://pfam.xfam.org/)

<!---
### Setup




```r
# Packages
library(MSnSet.utils)
library(dplyr)
library(kableExtra) # Table formatting

# Load the oca.set MSnSet
data("cptac_oca")

# DEA
res <- limma_a_b(oca.set,
                 model.str = "~ PLATINUM.STATUS", 
                 coef.str = "PLATINUM.STATUS")
```

There are no proteins with adjusted p-values less than 0.05, so we will use the unadjusted p-values to determine the differentially expressed proteins. (This is purely for the demonstration. You shouldn't actually do this.)

Of the 8103 proteins that were tested, [placeholder] are differentially expressed (based on unadjusted p-values). We need to convert the RefSeq protein identifiers to entrez gene identifiers if we want to perform ORA with the Reactome and KEGG databases. Some functions allow other keytypes such as REFSEQ, but we will use ENTREZID for everything to be consistent.


```r
# Map RefSeq protein IDs to Entrez Gene ID
conv_tbl <- MSnID::fetch_conversion_table("Homo sapiens",
                                          from = "REFSEQ",
                                          to = "ENTREZID")

res <- res %>% 
  tibble::rownames_to_column("REFSEQ") %>% 
  # Remove isoforms
  mutate(REFSEQ = gsub("(.*)\\..*", "\\1", REFSEQ)) %>%
  left_join(conv_tbl)
```

Now, we will create a vector for the genes of interest (the DEGs from `res`) and all genes that were tested. This is all we need for ORA.


```r
## Setup for ORA:
# The genes of interest are the unique entrez gene IDs with
# unadjusted p-values less than 0.05.
genes <- unique(res$ENTREZID[res$P.Value < 0.05]) # 359 Entrez Gene IDs

# The universe for the hypergeometric calculations.
# Sometimes called the vector of background genes.
universe <- unique(res$ENTREZID) # 7661 Entrez Gene IDs
```


### Gene Ontology

We will use `clusterProfiler::enrichGO` for ORA of the gene ontologies. For this example, we will use the biological process ontology.


```r
# Required packages
library(org.Hs.eg.db)
library(clusterProfiler)

ORA_BP <- enrichGO(gene = genes, 
                   universe = universe, 
                   OrgDb = "org.Hs.eg.db", 
                   keyType = "ENTREZID", 
                   ont = "BP",
                   pvalueCutoff = 0.05, 
                   qvalueCutoff = 1, 
                   pAdjustMethod = "BH", 
                   minGSSize = 1, 
                   readable = TRUE)

# Remove redundancy of enriched GO terms.
ORA_BP <- simplify(ORA_BP, cutoff = 0.85)

# Get the results data frame.
ORA_BP_df <- ORA_BP@result 
```


```r
# Format data frame for kable
ORA_BP_df$geneID <- strsplit(ORA_BP_df$geneID, split = "/")

# Keep first two genes from each geneID vector and
# replace the rest with ellipses.
ORA_BP_df$geneID <- unlist(lapply(ORA_BP_df$geneID, function(x)
  paste0(paste(head(x, 2), collapse = "/"), "...")
))

ORA_BP_df <- ORA_BP_df %>% 
  mutate(p.adjust = signif(p.adjust, 2),
         qvalue = signif(qvalue, 2))

# Top 4 enriched GO terms
kable(head(ORA_BP_df, 4), row.names = FALSE,
      caption = "Results from clusterProfiler::enrichGO") %>% 
  kable_styling(full_width = FALSE)
```

We can visualize the ORA results with a dotplot.


```r
# Dotplot of ORA results
enrichplot::dotplot(ORA_BP, showCategory = 10)
```

MSnSet.utils also has its own version of `dotplot` called `plot_enrichment`, though it is currently limited to objects of class `enrichResult`. One major benefit over `enrichplot::dotplot` is that the transformation for the color gradient (either "identity" or "log10") is automatically chosen based on the data. It also cuts back on the number of plot aesthetics (e.g. point size is not used) and uses sentence case for terms.


```r
plot_enrichment(ORA_BP, num.categories = 10, 
                labels = scales::label_scientific())
```


### Reactome


```r
library(ReactomePA)

ORA_reactome <- ReactomePA::enrichPathway(gene = genes,
                                          organism = "human", 
                                          pvalueCutoff = 0.05, 
                                          pAdjustMethod = "BH", 
                                          qvalueCutoff = 1, 
                                          universe = universe, 
                                          minGSSize = 1, 
                                          readable = TRUE)

ORA_reactome_df <- ORA_reactome@result
```


```r
# Format data frame for kable
ORA_reactome_df$geneID <- strsplit(ORA_reactome_df$geneID, split = "/")

# Keep first two genes from each geneID vector and
# replace the rest with ellipses.
ORA_reactome_df$geneID <- unlist(lapply(ORA_reactome_df$geneID, function(x)
  paste0(paste(head(x, 2), collapse = "/"), "...")
))

ORA_reactome_df <- ORA_reactome_df %>% 
  mutate(p.adjust = signif(p.adjust, 2),
         qvalue = signif(qvalue, 2))

# Top 4 enriched GO terms
kable(head(ORA_reactome_df, 4), row.names = FALSE,
      caption = "Results from ReactomePA::enrichPathway") %>%
  kable_styling(full_width = FALSE)
```


### KEGG


```r
ORA_KEGG <- clusterProfiler::enrichKEGG(gene = genes, 
                                        organism = "hsa", 
                                        keyType = "ncbi-geneid", 
                                        pvalueCutoff = 0.05, 
                                        qvalueCutoff = 1,
                                        pAdjustMethod = "BH",
                                        universe = universe, 
                                        minGSSize = 1, 
                                        use_internal_data = FALSE)

ORA_KEGG_db <- ORA_KEGG@result
```


```r
kable(head(ORA_KEGG_db, 4))
```



### Pfam
--->



