README
================

# ARscore

<!-- badges: start -->
<!-- badges: end -->

Phage immunoprecipitation sequencing (PhIP-Seq) is a massively
multiplexed method for quantifying antibody reactivity to libraries of
peptides. PhIP-seq analyses begin by identifying enriched antibody
reactivity to individual peptides. However, studies frequently require
understanding of aggregate reactivity to whole antigens or pathogens.
`ARscore` provides a standardized approach for calculating aggregate
reactivity to groups of peptides initially implemented with VirScan (a
PhIP-Seq library of viral peptides) to create a virus level reactivity
metric.

The ARscore method was applied to PhIP-Seq data where peptide
enrichments were determined with
[`edgeR`](https://bioconductor.org/packages/release/bioc/html/edgeR.html)’s
standard pipeline for identifying differential expression from read
count data[^1][^2][^3].

For more information, see the package vignette using
`browseVignettes("ARscore")`.

## Installation

### `ARscore`

To install the package:

``` r
if(!requireNamespace("remotes")) install.packages("remotes")

remotes::install_github("wmorgen1/ARscore")
```

To load the package:

``` r
library(ARscore)
```

[^1]: Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a
    Bioconductor package for differential expression analysis of digital
    gene expression data. Bioinformatics 26, 139-140

[^2]: McCarthy DJ, Chen Y and Smyth GK (2012). Differential expression
    analysis of multifactor RNA-Seq experiments with respect to
    biological variation. Nucleic Acids Research 40, 4288-4297

[^3]: Chen Y, Lun ATL, Smyth GK (2016). From reads to genes to pathways:
    differential expression analysis of RNA-Seq experiments using
    Rsubread and the edgeR quasi-likelihood pipeline. F1000Research 5,
    1438
