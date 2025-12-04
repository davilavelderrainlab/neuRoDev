--- 
title: "neuRoDev"
author: "Erik Bot & Asia Zonca"
date: "2025-12-04"
site: bookdown::bookdown_site
documentclass: book
bibliography: [book.bib, packages.bib]
# url: your book url like https://bookdown.org/yihui/bookdown
# cover-image: path to the social sharing image like images/cover.jpg
description: |
  This is the tutorial of the use of the neuRoDev package to explore the processes of corticogenesis.
biblio-style: apalike
csl: chicago-fullnote-bibliography.csl
---

# Summary
This is the tutorial of the use of the `neuRoDev` package from *doi* to explore the processes of corticogenesis. The following chapters allow the inspection of the resource compendium including both analyses shown in the article and additional examples and mode of use. We also included an interactive **eTrace** tool to investigate patterns of gene(s) expression instantaneously.

## Installation 
To install the `neuRoDev` package from *Github*:

``` r
install.packages("devtools")

devtools::install_github("https://github.com/davilavelderrainlab/neuRoDev")
```

The `neuRoDev` package uses SingleCellExperiment objects to store the corticogenesis, neurogenesis, and gliogenesis resources. 

The data used in this tutorial are available for download here.
