# Analysis Tools{#analysis}
## Loading necessary libraries and data

### Loading the right libraries




``` r
# functions package
library(neuRoDev)
# extra packages
library(ComplexHeatmap)
```

### Loading the networks necessary objects


``` r
corticogenesis_sce <- corticogenesis.sce(directory = '~/Downloads')
neurogenesis_sce <- neurogenesis.sce(directory = '~/Downloads')
gliogenesis_sce <- gliogenesis.sce(directory = '~/Downloads')
```

## eTraces

We developed an intuitive way to directly visualize a molecular phenotype score in each cluster of our reference networks, which we called **eTrace** (expression Trace). We provide here an [interactive tool](#interactive) to directly assess the (log-normalized) expression of input genes and gene sets, as well as the expression enrichment trends (more statistically robust). 

The function used to visualize eTraces is **plot_eTrace**. The function requires few inputs, with the only mandatory one being the reference network (`net`):

-   `net`= the reference network to use.
-   `genes`= a specific gene or set of genes from which to derive the score to plot. If more genes are given, the output will be an averaged score of those genes. It defaults to NULL, as the user can directly provide a score per cluster (see `score`). If no genes are given and a score matrix is selected, all genes in the score matrix will be considered.
-   `score`= it can be a vector, with one value per cluster, or a matrix. If a matrix is given together with genes (see `genes`), only the values of the selected genes that are present in the matrix will be shown. It defaults to the log-normalized expression values contained in `net` under `logcounts`.
-   `expression_enrichment`= it defines whether to compute expression enrichment of the given genes, a more statistically robust way of looking at the expression. If TRUE, the score that will be used will be the expression enrichment score, regardless on what is put in `score`. It defaults to FALSE.
-   other plotting and fine-tuning inputs, see **?plot_eTrace** for further information.


``` r
plot_eTrace(corticogenesis_sce, 
            genes = 'SST', 
            main = 'SST - no enrichment')
```

<div class="figure">
<img src="3-Analysis_Tools_files/figure-html/ch3-fig1-1.png" alt="Single gene eTrace in corticogenesis." width="90%" />
<p class="caption">(\#fig:ch3-fig1)Single gene eTrace in corticogenesis.</p>
</div>


``` r
plot_eTrace(corticogenesis_sce, 
            genes = 'PAX6', 
            expression_enrichment = TRUE, 
            main = 'PAX6 - enrichment')
```

<div class="figure">
<img src="3-Analysis_Tools_files/figure-html/ch3-fig2-1.png" alt="Single gene expression enrichment eTrace in corticogenesis." width="90%" />
<p class="caption">(\#fig:ch3-fig2)Single gene expression enrichment eTrace in corticogenesis.</p>
</div>

It is possible to inspect not only single genes but also gene sets. For example, we can visualize the expression of curated preferentially expressed genes in the reference subclasses. The objects are available for download here. 


``` r
corticogenesis_pe_genes <- readRDS("~/Downloads/corticogenesis_subclass_preferential_genes.rds")
neurogenesis_pe_genes <- readRDS("~/Downloads/neurogenesis_subclass_preferential_genes.rds")
gliogenesis_pe_genes <- readRDS("~/Downloads/gliogenesis_subclass_preferential_genes.rds")
```


``` r
plot_eTrace(corticogenesis_sce, 
            genes = corticogenesis_pe_genes$Oligo, 
            main = 'Oligodendrocyte preferentially expressed genes')
```

<div class="figure">
<img src="3-Analysis_Tools_files/figure-html/ch3-fig3-1.png" alt="Gene set eTrace in corticogenesis." width="90%" />
<p class="caption">(\#fig:ch3-fig3)Gene set eTrace in corticogenesis.</p>
</div>


``` r
plot_eTrace(corticogenesis_sce, 
            genes = corticogenesis_pe_genes$Opc, 
            expression_enrichment = TRUE, 
            main = 'OPC preferentially expressed genes')
```

<div class="figure">
<img src="3-Analysis_Tools_files/figure-html/ch3-fig4-1.png" alt="Genes set expression enrichment eTrace in corticogenesis." width="90%" />
<p class="caption">(\#fig:ch3-fig4)Genes set expression enrichment eTrace in corticogenesis.</p>
</div>

Any kind of score can be visualized with the eTrace. As a further example, we can visualize the preferential expression of Gene Ontology genesets. 

We have already computed preferential expression of Gene Ontology Biological Processes (BP), Molecular Functions (MF), and Cellular Components (CC), which can be downloaded here. 
Each object is a list containing preferential expression scores in one of the three reference networks in the three ontologies (BP, MF, CC). Each element of each list contains the activity (`activity`) derived from Gene Set Variation Analysis (one value per gene set in each cluster) and the preferential expression scores (`preferential`; one value per gene set in each subclass).


``` r
corticogenesis_preferential_GO <- readRDS('~/Downloads/corticogenesis_preferential_GO.rds')
neurogenesis_preferential_GO <- readRDS('~/Downloads/neurogenesis_preferential_GO.rds')
gliogenesis_preferential_GO <- readRDS('~/Downloads/gliogenesis_preferential_GO.rds')
```


``` r
plot_eTrace(corticogenesis_sce, 
            score = corticogenesis_preferential_GO$GO_Biological_Process_2025$activity["Glycogen Biosynthetic Process (GO:0005978)",], 
            main = 'Glycogen Biosynthesis')
```

<div class="figure">
<img src="3-Analysis_Tools_files/figure-html/ch3-fig5-1.png" alt="Gene set activity score eTrace in corticogenesis." width="90%" />
<p class="caption">(\#fig:ch3-fig5)Gene set activity score eTrace in corticogenesis.</p>
</div>

The expression enrichment values can be obtained by using the function **get_eTrace**, which requires as input: 
  - `net`= the reference network.
  - `genes`= the genes to use.
  - `nRand`= number of random genes to use as a comparison. Defaults to 100.

Additionally, it is possible to look at within-lineage specific patterns of expression by subsetting the networks. For example, here we show the expression patterns of genes inside the astrogenesis- and oligodendrogenesis-specific trajectories of differentiation:

``` r
astro_sce <- gliogenesis_sce[,-c(grep(gliogenesis_sce$SubClass, pattern="OPC"),
  grep(gliogenesis_sce$SubClass, pattern="Oli"))] #removes everything that is oligodendroglia-related

oli_sce <- gliogenesis_sce[,c(grep(gliogenesis_sce$SubClass, pattern="OPC"),
  grep(gliogenesis_sce$SubClass, pattern="Oli"))] #keeps only oligodendroglia-related
```


``` r
plot_eTrace(astro_sce,
            genes = "SPARCL1",
            expression_enrichment = T,
            main = "SPARCL1 in astrogenesis")

plot_eTrace(oli_sce,
            genes = "MAG",
            expression_enrichment = T,
            main = "MAG in oligodendrogenesis")
```



<div class="figure">
<img src="3-Analysis_Tools_files/figure-html/ch3-fig6-1.png" alt="Single gene expression enrichment eTrace within astrogenesis." width="90%" />
<p class="caption">(\#fig:ch3-fig6)Single gene expression enrichment eTrace within astrogenesis.</p>
</div>

<div class="figure">
<img src="3-Analysis_Tools_files/figure-html/ch3-fig7-1.png" alt="Single gene expression enrichment eTrace within oligodendrogenesis." width="90%" />
<p class="caption">(\#fig:ch3-fig7)Single gene expression enrichment eTrace within oligodendrogenesis.</p>
</div>


## Interactive eTrace{#interactive}
To explore the normalized expression and expression enrichment levels of (single) genes it is possible to use this interactive tool and observe in a glance patterns over ~70 years of life. 
Below you can access with one click to the three different resources and give as input a single gene or a gene set and look at the patterns of expression in time and through subclasses.  

<a href="https://erikbot.shinyapps.io/etraceshinycortico/" target="_blank">
  🔍 Click here to interactively visualize the corticogenesis eTrace
</a>

<a href="https://erikbot.shinyapps.io/etraceshinyneuro/" target="_blank">
  🔍 Click here to interactively visualize the neurogenesis eTrace
</a>

<a href="https://erikbot.shinyapps.io/etraceshinyglio/" target="_blank">
  🔍 Click here to interactively visualize the gliogenesis eTrace
</a>


## Expression enrichment across stage and subclass

An alternative and compact visualization of expression enrichment in both subclasses and stages can be obtained by using the function **plot_eMatrix**. The inputs are listed below, and are the same as those needed by the **get_eTrace** function:

  - `net`= the reference network to use.
  - `genes`= a specific gene or set of genes from which to derive the expression enrichment.
  - `nRand`= the number of random sets to use for the expression enrichment calculation. Defaults to 100.

This returns an **eMatrix** showing on the columns subclasses, on the rows stages, and the values represent enrichment scores. 


``` r
plot_eMatrix(net = corticogenesis_sce, 
             genes = 'VIP')
```

<div class="figure" style="text-align: center">
<img src="3-Analysis_Tools_files/figure-html/ch3-fig8-1.png" alt="VIP eMatrix in corticogenesis." width="441.6" />
<p class="caption">(\#fig:ch3-fig8)VIP eMatrix in corticogenesis.</p>
</div>

<details>
<summary><strong>Show the code (plotting eMatrix)</strong></summary>

``` r
plot_eMatrix(net = neurogenesis_sce, 
             genes = 'NEUROD6')
plot_eMatrix(net = gliogenesis_sce, 
             genes = 'GFAP')
plot_eMatrix(net = corticogenesis_sce, 
             genes = corticogenesis_pe_genes$L6NIT)
```
</details>

<div class="figure" style="text-align: center">
<img src="3-Analysis_Tools_files/figure-html/ch3-fig9-1.png" alt="NEUROD6 eMatrix in neurogenesis." width="528" />
<p class="caption">(\#fig:ch3-fig9)NEUROD6 eMatrix in neurogenesis.</p>
</div>

<div class="figure" style="text-align: center">
<img src="3-Analysis_Tools_files/figure-html/ch3-fig10-1.png" alt="GFAP eMatrix in gliogenesis." width="432" />
<p class="caption">(\#fig:ch3-fig10)GFAP eMatrix in gliogenesis.</p>
</div>

<div class="figure" style="text-align: center">
<img src="3-Analysis_Tools_files/figure-html/ch3-fig11-1.png" alt="L6NIT preferentially expressed genes eMatrix in corticogenesis." width="528" />
<p class="caption">(\#fig:ch3-fig11)L6NIT preferentially expressed genes eMatrix in corticogenesis.</p>
</div>

It is possible to obtain the `eMatrix` itself by using the function **get_eMatrix**, which uses the same inputs as **plot_eMatrix**.


## Visualization of enrichment in the network

As described in [Chapter Network exploration](#network), we can visualize cluster-wise scores directly in the UMAP at the network level. This can be done with the **plotNetworkScore** function, setting `expression_enrichment` to TRUE.


``` r
plotNetworkScore(net = corticogenesis_sce, 
                 genes = corticogenesis_pe_genes$Oligo, 
                 expression_enrichment = TRUE,
                 main = "Enrichment: Oligo preferential genes")
```

<div class="figure" style="text-align: center">
<img src="3-Analysis_Tools_files/figure-html/ch3-fig12-1.png" alt="Gene expression enrichment in the corticogenesis network." width="672" />
<p class="caption">(\#fig:ch3-fig12)Gene expression enrichment in the corticogenesis network.</p>
</div>

## Preferential expression
The reference networks also contain the pseudobulks of each subclass and each stage, along with the preferential expression profiles computed between subclasses/stages. We can inspect the genes with the highest preferential expression scores in each subclass. 

<details>
<summary><strong>Show the code (plotting heatmaps)</strong></summary>

``` r
top_genes_cortico <- unique(rownames(corticogenesis_sce)[apply(assays(corticogenesis_sce@metadata$subclass_psb)[['preferential']], 2, function(i) {
  order(i, decreasing =  T)[1:5]})])

h_cortico <- Heatmap(t(assays(corticogenesis_sce@metadata$subclass_psb)[['preferential']][top_genes_cortico, levels(corticogenesis_sce$SubClass)]), 
  name = 'preferential\nexpression', 
  cluster_columns = T,
  cluster_rows = F,
  height = grid::unit(ncol(assays(corticogenesis_sce@metadata$subclass_psb)[['preferential']])*4, 'mm'), 
  width = grid::unit(length(top_genes_cortico)*4, 'mm'),
  rect_gp = gpar(color = 'black', lwd = 0.5))

draw(h_cortico, , heatmap_legend_side = 'left')
```


``` r
top_genes_neuro <- unique(rownames(neurogenesis_sce)[apply(assays(neurogenesis_sce@metadata$subclass_psb)[['preferential']], 2, function(i) {
  order(i, decreasing =  T)[1:5]})])

h_neuro <- Heatmap(t(assays(neurogenesis_sce@metadata$subclass_psb)[['preferential']][top_genes_neuro, levels(neurogenesis_sce$SubClass)]), 
  name = 'preferential expression', 
  cluster_columns = T,
  cluster_rows = F,
  height = grid::unit(ncol(assays(neurogenesis_sce@metadata$subclass_psb)[['preferential']])*4, 'mm'), 
  width = grid::unit(length(top_genes_neuro)*4, 'mm'), 
  rect_gp = gpar(color = 'black', lwd = 0.5))

draw(h_neuro, , heatmap_legend_side = 'left')
```


``` r
top_genes_glio <- unique(rownames(gliogenesis_sce)[apply(assays(gliogenesis_sce@metadata$subclass_psb)[['preferential']], 2, function(i) {
  order(i, decreasing =  T)[1:5]})])

h_glio <- Heatmap(t(assays(gliogenesis_sce@metadata$subclass_psb)[['preferential']][top_genes_glio, levels(gliogenesis_sce$SubClass)]), 
  name = 'preferential expression', 
  cluster_columns = T,
  cluster_rows = F,
  height = grid::unit(ncol(assays(gliogenesis_sce@metadata$subclass_psb)[['preferential']])*4, 'mm'),
  width = grid::unit(length(top_genes_glio)*4, 'mm'), 
  rect_gp = gpar(color = 'black', lwd = 0.5))

draw(h_glio, , heatmap_legend_side = 'left')
```
</details>

<div class="figure">
<img src="3-Analysis_Tools_files/figure-html/ch3-fig13-1.png" alt="Top preferential genes in corticogenesis." width="92%" />
<p class="caption">(\#fig:ch3-fig13)Top preferential genes in corticogenesis.</p>
</div>

<div class="figure">
<img src="3-Analysis_Tools_files/figure-html/ch3-fig14-1.png" alt="Top preferential genes in neurogenesis." width="92%" />
<p class="caption">(\#fig:ch3-fig14)Top preferential genes in neurogenesis.</p>
</div>

<div class="figure">
<img src="3-Analysis_Tools_files/figure-html/ch3-fig15-1.png" alt="Top preferential genes in gliogenesis." width="95%" />
<p class="caption">(\#fig:ch3-fig15)Top preferential genes in gliogenesis.</p>
</div>

At the same time, we computed the preferential expression scores for Gene Ontology gene sets. As an example, we can visualize the preferential expression scores of Gene Ontology Biological Processes.
<details>
<summary><strong>Show the code (plotting heatmaps)</strong></summary>

``` r
top_GO_bp_cortico_list <- apply(corticogenesis_preferential_GO$GO_Biological_Process_2025$preferential, 2, function(i) {
  rownames(corticogenesis_preferential_GO$GO_Biological_Process_2025$preferential)[order(i, decreasing = TRUE)[seq(1,2)]]
  }
) #selecting top 2 processes by subclass

top_GO_bp_cortico <- corticogenesis_preferential_GO$GO_Biological_Process_2025$preferential[unique(as.vector(top_GO_bp_cortico_list)),]

rownames(top_GO_bp_cortico) <- sub("\\s*\\(GO:\\d+\\)$", "", rownames(top_GO_bp_cortico))

h <- Heatmap(top_GO_bp_cortico[,levels(corticogenesis_sce$SubClass)], 
             name = 'preferential\nexpression', 
             width = grid::unit(ncol(top_GO_bp_cortico)*4, 'mm'), 
             height = grid::unit(nrow(top_GO_bp_cortico)*4, 'mm'), 
             cluster_columns = FALSE, 
             rect_gp = gpar(color = 'black', lwd = 0.5))
draw(h, heatmap_legend_side = 'left')
```


``` r
top_GO_bp_neuro_list <- apply(neurogenesis_preferential_GO$GO_Biological_Process_2025$preferential, 2, function(i) {
  rownames(neurogenesis_preferential_GO$GO_Biological_Process_2025$preferential)[order(i, decreasing = TRUE)[seq(1,2)]]
  }
)

top_GO_bp_neuro <- neurogenesis_preferential_GO$GO_Biological_Process_2025$preferential[unique(as.vector(top_GO_bp_neuro_list)),]

rownames(top_GO_bp_neuro) <- sub("\\s*\\(GO:\\d+\\)$", "", rownames(top_GO_bp_neuro))

h_neuro <- Heatmap(top_GO_bp_neuro[,levels(neurogenesis_sce$SubClass)], 
                   name = 'preferential\nexpression', 
                   width = grid::unit(ncol(top_GO_bp_neuro)*4, 'mm'), 
                   height = grid::unit(nrow(top_GO_bp_neuro)*4, 'mm'), 
                   cluster_columns = FALSE, 
                   rect_gp = gpar(color = 'black', lwd = 0.5))
draw(h_neuro, heatmap_legend_side = 'left')
```


``` r
top_GO_bp_glio_list <- apply(gliogenesis_preferential_GO$GO_Biological_Process_2025$preferential, 2, function(i) {
  rownames(gliogenesis_preferential_GO$GO_Biological_Process_2025$preferential)[order(i, decreasing = TRUE)[seq(1,2)]]
  }
)

top_GO_bp_glio <- gliogenesis_preferential_GO$GO_Biological_Process_2025$preferential[unique(as.vector(top_GO_bp_glio_list)),]

rownames(top_GO_bp_glio) <- sub("\\s*\\(GO:\\d+\\)$", "", rownames(top_GO_bp_glio))

h_glio <- Heatmap(top_GO_bp_glio[,levels(gliogenesis_sce$SubClass)], 
                  name = 'preferential\nexpression', 
                  width = grid::unit(ncol(top_GO_bp_glio)*4, 'mm'), 
                  height = grid::unit(nrow(top_GO_bp_glio)*4, 'mm'), 
                  cluster_columns = FALSE, 
                  rect_gp = gpar(color = 'black', lwd = 0.5))
draw(h_glio, heatmap_legend_side = 'left')
```
</details>




<div class="figure">
<img src="3-Analysis_Tools_files/figure-html/ch3-fig16-1.png" alt="Top gene ontology biological processes in corticogenesis." width="90%" />
<p class="caption">(\#fig:ch3-fig16)Top gene ontology biological processes in corticogenesis.</p>
</div>

<div class="figure">
<img src="3-Analysis_Tools_files/figure-html/ch3-fig17-1.png" alt="Top gene ontology biological processes in neurogenesis." width="90%" />
<p class="caption">(\#fig:ch3-fig17)Top gene ontology biological processes in neurogenesis.</p>
</div>

<div class="figure">
<img src="3-Analysis_Tools_files/figure-html/ch3-fig18-1.png" alt="Top gene ontology biological processes in gliogenesis." width="90%" />
<p class="caption">(\#fig:ch3-fig18)Top gene ontology biological processes in gliogenesis.</p>
</div>
To focus on corticogenesis-relevant processes we have also manually curated a list of Gene Ontology Biological Processes for the corticogenesis and neurogenesis networks, available for download here.
