neuRoDev is an R package that provides analysis tools to make use of the 
integrated transcriptomic references of cortex development, described in 
[Zonca&, Bot&, and Davila-Velderrain, 2025].

All the datasets needed as input are stored in the [Figshare database].
Functions in the package can install the datasets automatically (see
[Network exploration]).

# Analysis tools

neuRoDev contains functions to perform:

- Reference networks exploration ([Network exploration])
- Expression enrichment and specificity ([Analysis tools])
- Single-cell RNAseq datasets mapping ([Mapping scRNAseq])
- Bulk RNAseq datasets mapping ([Mapping bulkRNAseq])

In the [Tutorial], we provide an extensive description and the code to perform
all the mentioned analyses.

# Basic installation

To install `neuRoDev` from GitHub:

```{r}
install.packages("devtools")

devtools::install_github("https://github.com/davilavelderrainlab/neuRoDev")
```

# Bug report

Please use the [issues] to submit bug reports.

# Reference

If you use `neuRoDev` in your work, please cite

> **Developmental single-cell transcriptomic networks dissect neuronal and glial genesis and maturation in the human cortex**
>
> Asia Zonca&, Erik Bot& & José Davila-Velderrain
>
> _Journal_ Date. doi: [doi](https://github.com/davilavelderrainlab/neuRoDev).

[Zonca&, Bot&, and Davila-Velderrain, 2025]: https://github.com/davilavelderrainlab/neuRoDev
[Figshare database]: https://github.com/davilavelderrainlab/neuRoDev
[Network exploration]: https://github.com/davilavelderrainlab/neuRoDev
[Analysis tools]: https://github.com/davilavelderrainlab/neuRoDev
[Gene set networks]: https://github.com/davilavelderrainlab/neuRoDev
[Mapping scRNAseq]: https://github.com/davilavelderrainlab/neuRoDev
[Mapping bulkRNAseq]: https://github.com/davilavelderrainlab/neuRoDev
[Tutorial]: https://github.com/davilavelderrainlab/neuRoDev
[article]: https://github.com/davilavelderrainlab/neuRoDev
[issues]: https://github.com/davilavelderrainlab/neuRoDev/issues
