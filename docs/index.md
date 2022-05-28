--- 
title: "Proteomics Data Analysis in R/Bioconductor"
author: "Tyler Sagendorf"
date: "May 27, 2022"
site: bookdown::bookdown_site
output: bookdown::gitbook
documentclass: book
bibliography:
- references.bib
- r-packages.bib
csl: begell-house-chicago-author-date.csl
link-citations: yes
nocite: '@*'
description: This is a tutorial for proteomics data analysis in R that utilizes packages
  developed by researchers at PNNL and from Bioconductor.
---

# Welcome! {-}

This tutorial is very much a work-in progress. Even sections that appear finished are likely to be changed. I will update this when significant progress is made. Thank you for your patience.




It is highly recommended to review the resources below before continuing with the rest of the tutorial.


* Proteomics Overview
    * <a href="https://pubs.acs.org/doi/10.1021/cr3003533">Protein Analysis by Shotgun/Bottom-up Proteomics</a>
    * <a href="https://link.springer.com/book/10.1007%2F978-3-319-41448-5">Modern Proteomics – Sample Preparation, Analysis and Practical Applications</a>
    * <a href="https://dx.doi.org/10.1214%2F10-AOAS341">Liquid Chromatography Mass Spectrometry-Based Proteomics: Biological and Technological Aspects</a>

<!---
* High-Performance Liquid Chromatography (HPLC)
--->

* Mass Spectrometry
    * <a href="https://warwick.ac.uk/fac/sci/lifesci/research/sigtraf/animations/">Warwick School of Life Sciences Teaching Animations</a>
    * <a href="https://doi.org/10.2144/05384te01">Tandem Mass Spectrometry for Peptide and Protein Sequence Analysis</a>
    * <a href="https://escholarship.org/uc/item/8tt6h3jt">Maestro: Comprehensive, Multi-Stage Spectrum Identification in Protein Mass Spectrometry</a>
    * <a href="https://youtu.be/Esf1EqzyQZc">Searching databases for protein identification - part 1</a> (YouTube video)
    * <a href="https://www.youtube.com/watch?v=v8EsEWwrJWs">Mass spectrometry for proteomics - part one</a> (YouTube video)
    * <a href="http://www.ncbi.nlm.nih.gov/pmc/articles/pmc1853331/">Electrospray Ionisation Mass Spectrometry: Principles and Clinical Applications</a>


* <a href="https://dms2.pnl.gov/">PNNL's Data Management System (DMS)</a>


* <a href="https://prismwiki.pnl.gov/wiki/Category:Glossary" title="PRISMWiki Glossary">Integrative Omics PRISMWiki</a>


* Universal Protein Resource (<a href="https://www.uniprot.org/help/about" title="About UniProt">UniProt</a>): protein sequence and annotation data


* False Discovery Rate (FDR)
    * <a href="https://dx.doi.org/10.1002%2Fpmic.201500431">How to talk about protein‐level false discovery rates in shotgun proteomics</a>
    * <a href="https://pubs.acs.org/doi/10.1021/pr700739d">Posterior Error Probabilities and False Discovery Rates: Two Sides of the Same Coin</a>
    * <a href="https://www.bioinfor.com/fdr-tutorial/">False Discovery Rate: PEAKS FDR Estimation</a>
    * <a href="https://dx.doi.org/10.1186%2F1471-2105-13-S16-S2">False discovery rates in spectral identification</a>


* <a href = "https://www.rstudio.com/resources/cheatsheets/">RStudio Cheatsheets</a>


* Pattern matching with regular expressions
    * <a href = "https://r4ds.had.co.nz/strings.html">R for Data Science: Strings</a>
    * <a href="https://regexone.com/">RegexOne: Learn Regular Expressions with simple, interactive exercises.</a>


