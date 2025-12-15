# neuRoDev
neuRoDev is a computational resource with integrative networks, transcriptomes, and analytical tools to study the development of neuronal and glial cells in human

neuRoDev is decribed in the paper [Zonca, Bot, and Davila-Velderrain, 2025].

This repository contains the neuRoDev R package

Accompanying data resources can be found [here]

Follow installation functions in the package to automatically download and use all data sources (see Network exploration).

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
> _Journal_ Date. doi: [doi](https://github.com/davilavelderrainlab/neuRoDev).

[Zonca, Bot, and Davila-Velderrain, 2025]: https://github.com/davilavelderrainlab/neuRoDev
[here]: https://github.com/davilavelderrainlab/neuRoDev
[Network exploration]: https://github.com/davilavelderrainlab/neuRoDev
[Analysis tools]: https://github.com/davilavelderrainlab/neuRoDev
[Mapping scRNAseq]: https://github.com/davilavelderrainlab/neuRoDev
[Mapping bulkRNAseq]: https://github.com/davilavelderrainlab/neuRoDev
[corticogenesis]: https://erikbot.shinyapps.io/etraceshinycortico/
[neurogenesis]: https://erikbot.shinyapps.io/etraceshinyneuro/
[gliogenesis]: https://erikbot.shinyapps.io/etraceshinyglio/
[Tutorial]: https://github.com/davilavelderrainlab/neuRoDev
[article]: https://github.com/davilavelderrainlab/neuRoDev
[issues]: https://github.com/davilavelderrainlab/neuRoDev/issues
