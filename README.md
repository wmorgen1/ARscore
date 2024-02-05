
# ARscore

<!-- badges: start -->
<!-- badges: end -->

Phage ImmunoPrecipitation Sequencing (PhIP-Seq) is a massively
multiplexed method for quantifying antibody reactivity to libraries of
peptides. PhIP-seq analyses begin by identifying enriched antibody
reactivity to individual peptides. However, studies frequently require
understanding of aggregate reactivity to whole antigens or pathogens.
`ARscore` provides a standardized approach for calculating aggregate
reactivity to groups of peptides.

`ARscore` generates aggregate reactivity scores (ARscores) by comparing
the average fold change of a group of peptides to distributions of
average fold change from randomly selected peptides, using
[`fitdistrplus`](https://cran.r-project.org/web/packages/fitdistrplus/index.html)[^1]
and
[`limma`](https://bioconductor.org/packages/release/bioc/html/limma.html)[^2].
To expedite computation and evaluate extreme reactivity scores, random
distributions are modeled as gamma distributions with paramaters that
change regularly with the number of randomly selected peptides.

ARscore was initially implemented with VirScan (a PhIP-Seq library of
viral peptides) to create a virus level reactivity metric (Viral
ARscore, VARscore), using peptide enrichments determined with
[`edgeR`](https://bioconductor.org/packages/release/bioc/html/edgeR.html)’s
standard pipeline for identifying differential expression from read
count data[^3][^4][^5].

Input PhIP-Seq data requires Larman Lab naming conventions for mock IP
controls, samples, and peptide annotations. Peptide grouping is
currently based on the taxon_species annotation column. Custom peptide
groupings can be achieved by replacing this column.

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

[^1]: Delignette-Muller ML, Dutang C (2015). fitdistrplus: An R Package
    for Fitting Distributions. Journal of Statistical Software, 64(4),
    1–34.

[^2]: Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015).
    limma powers differential expression analyses for RNA-sequencing and
    microarray studies. Nucleic Acids Research, 43(7), e47.

[^3]: Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a
    Bioconductor package for differential expression analysis of digital
    gene expression data. Bioinformatics 26, 139-140.

[^4]: McCarthy DJ, Chen Y and Smyth GK (2012). Differential expression
    analysis of multifactor RNA-Seq experiments with respect to
    biological variation. Nucleic Acids Research 40, 4288-4297.

[^5]: Chen Y, Lun ATL, Smyth GK (2016). From reads to genes to pathways:
    differential expression analysis of RNA-Seq experiments using
    Rsubread and the edgeR quasi-likelihood pipeline. F1000Research 5,
    1438.
