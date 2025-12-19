# neuRoDev
neuRoDev is a computational resource with integrative networks, transcriptomes, and analytical tools to study the development of neuronal and glial cells in human

neuRoDev is decribed in the paper [Zonca, Bot, and Davila-Velderrain, 2025].

This repository contains the neuRoDev R package

Accompanying data resources can be found [here]

Follow installation functions in the package to automatically download and use all data sources (see [Tutorial]).

# Overview

in neuRoDev you can:

  - explore integrative neuronal and glial summary networks ([Network exploration])
  - analyze temporal and cellular variability with eTrace analysis ([Analysis tools])
  - map and interpret query transcriptomic data ([Mapping scRNAseq]) ([Mapping bulkRNAseq])

follow neuRoDev's [Tutorial] for extensive explanations and code examples

# Basic installation

To install `neuRoDev` from GitHub:

```{r}
install.packages("devtools")

devtools::install_github("https://github.com/davilavelderrainlab/neuRoDev", dependencies = TRUE)
```

# Bug report

Please use the [issues] to submit bug reports.

# Reference

If you use `neuRoDev` in your work, please cite

> **NeuRoDev resolves lifelong temporal and cellular variation in human cortical gene expression**
>
> Asia Zonca, Erik Bot & José Davila-Velderrain
>
> _bioRxiv_ 17/12/2025. doi: [https://doi.org/10.64898/2025.12.16.694589](https://doi.org/10.64898/2025.12.16.694589).

[Zonca, Bot, and Davila-Velderrain, 2025]: https://www.biorxiv.org/content/10.64898/2025.12.16.694589v1
[here]: https://doi.org/10.6084/m9.figshare.30885428
[Network exploration]: https://davilavelderrainlab.github.io/neuRoDev/network.html
[Analysis tools]: https://davilavelderrainlab.github.io/neuRoDev/analysis.html
[Mapping scRNAseq]: https://davilavelderrainlab.github.io/neuRoDev/mapping-sc.html
[Mapping bulkRNAseq]: https://davilavelderrainlab.github.io/neuRoDev/mapping-bulk.html
[Tutorial]: https:///davilavelderrainlab.github.io/neuRoDev
[issues]: https://github.com/davilavelderrainlab/neuRoDev/issues
