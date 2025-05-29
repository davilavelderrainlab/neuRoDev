#' Create the gene set network
#'
#' @param geneset_list The list of annotated gene sets
#' @param pseudo The pseudobulk in form of a SingleCellExperiment object (it
#' needs to have the logcounts attribute)
#' @param net_df The reference network dataframe (the one which contains
#' the correlation values with the reference signatures, obtained from
#' `reference_signatures_correlation`)
#' @param stages The attribute in which to find the time information (stages) in
#' the `net_df`. Defaults to 'Stages'.
#' @param resolution The resolution for clustering in igraph::cluster_leiden.
#' Defaults to 4.
#' @param padj_threshold The adjusted pvalue threshold to consider one
#' geneset to be significant. Defaults to 0.05
#' @param gsea If the GSEA was already performed, it can be given here and not
#' recomputed.
#' @param between_class_pseudo_enrichment If an additional pseudobulk enrichment
#' is performed, for example between classes, it can be given here as a single
#' vector for the wanted class.
#' @param between_class_gsea If the GSEA on the additional pseudobulk enrichment
#' was already performed before, it can be given here and not recomputed.
#' @param pseudo_de If the pseudobulk enrichment between time points (stages)
#' was already computed, it can be given here and not recomputed.
#' @param max.words The maximum number of keywords per cluster to be returned.
#' @param stages_palette The palette for the visualization of the stages.
#' @param ... Additional `fgsea` parameters
#'
#' @return A list that contains the GSEA(s), the geneset network and an
#' geneset dataframe, with the genesets, geneset groups, the pvalue
#' adjusted of the GSEA(s) (one column per GSEA), the stage in which the group
#' has the maximum expression, the leading edges of each geneset and the top
#' 10 most frequent genes for each geneset group.
#' @export
#'
#' @examples
#' geneset_list <- list('A' = c('Gene-1', 'Gene-2', 'Gene-3'),
#' 'B' = c('Gene-1', 'Gene-2', 'Gene-4'),
#' 'C' = c('Gene-3', 'Gene-4', 'Gene-1'))
#' pseudo <- SingleCellExperiment::SingleCellExperiment(
#' assays = list(logcounts = matrix(sample(seq(1,100),
#' 4*10,
#' replace = TRUE),
#' nrow = 4)))
#' rownames(pseudo) <- paste0('Gene-', seq(1, nrow(pseudo)))
#' colnames(pseudo) <- paste0('Column-', seq(1, ncol(pseudo)))
#' S <- FC_signatures(matrix(runif(200,0,10), ncol = 10))
#' rownames(S) <- paste0('Gene-', seq(1, dim(S)[1]))
#' refS <- FC_signatures(matrix(runif(200,0.1,7), ncol = 10))
#' colnames(refS) <- paste0('Reference-', seq(1, dim(refS)[2]))
#' rownames(refS) <- paste0('Gene-', seq(1, dim(refS)[1]))
#' net_df <- reference_signatures_correlation(S, refS)
#' net_df$Stages <- paste0('Stage.', c(1, 1, 2, 2, 2))
#' stages_palette <- c('Stage.1' = 'red', 'Stage.2' = 'blue')
#' get_geneset_network(geneset_list, pseudo, net_df)
get_geneset_network <- function(geneset_list,
                                pseudo,
                                net_df,
                                stages_palette,
                                stages = 'Stages',
                                resolution = 4,
                                padj_threshold = 0.05,
                                gsea = NULL,
                                between_class_pseudo_enrichment = NULL,
                                between_class_gsea = NULL,
                                pseudo_de = NULL,
                                max.words = 10,
                                ...) {

  lookup <- stats::setNames(rep(TRUE, length(rownames(pseudo))), rownames(pseudo))

  fast_filter <- function(vec, lookup) {
    vec[!is.na(lookup[as.character(vec)])]
  }

  geneset_list <- lapply(geneset_list, fast_filter, lookup = lookup)
  geneset_list <- geneset_list[which(unlist(lapply(geneset_list, length)) > 2)]

  net_df <- net_df[which(net_df$Cluster %in% colnames(pseudo)),]

  pseudo <- pseudo[,match(net_df$Cluster, colnames(pseudo))]

  if(is.vector(stages) & length(stages) == 1) {
    stages <- net_df[,stages]
  }

  if(is.null(pseudo_de)) {
    pseudo_de <- get_pseudobulk_preferential_expression(SingleCellExperiment::logcounts(pseudo),
                                                        groups = make.names(stages))
  }

  if(!is.null(gsea)) {
    fgseas_variable <- gsea
  } else {
    fgseas_variable <- apply(pseudo_de, 2, function(v) {

      names(v) <- rownames(pseudo_de)

      p <- fgsea::fgsea(geneset_list, stats = v, ...)

      p <- p[order(p$padj),]

      return(p)

    })
  }

  fgseas_variable_pathways <- lapply(fgseas_variable, function(i) {
    return(i$pathway[which(i$padj < padj_threshold & i$NES > 0)])
  })

  if(length(unlist(fgseas_variable_pathways)) == 0) {
    return('No significant pathways found in any comparison')
  }

  fgseas_variable_leadingedges <- lapply(fgseas_variable, function(i) {
    le <- i$leadingEdge[which(i$padj < padj_threshold & i$NES > 0)]
    names(le) <- i$pathway[which(i$padj < padj_threshold & i$NES > 0)]
    return(le)
  })

  fgseas_variable_leadingedges_no_names <- fgseas_variable_leadingedges
  names(fgseas_variable_leadingedges_no_names) <- NULL

  all_pathways <- unlist(fgseas_variable_pathways)
  GeneList <- unlist(fgseas_variable_leadingedges_no_names, recursive = FALSE)

  if(is.null(between_class_gsea) & !is.null(between_class_pseudo_enrichment)) {

    between_class_gsea <- fgsea::fgsea(geneset_list,
                                       stats = between_class_pseudo_enrichment,
                                       ...)

    between_class_gsea <- between_class_gsea[order(between_class_gsea$padj),]

  }

  if(!is.null(between_class_gsea)) {

    fgseas_invariable_pathways <-  between_class_gsea$pathway[which(between_class_gsea$padj < padj_threshold & between_class_gsea$NES > 0)]
    all_pathways <- c(fgseas_invariable_pathways, all_pathways)

    fgseas_invariable_leadingedges <- between_class_gsea$leadingEdge[which(between_class_gsea$padj < padj_threshold & between_class_gsea$NES > 0)]
    names(fgseas_invariable_leadingedges) <- fgseas_invariable_pathways

    GeneList <- c(fgseas_invariable_leadingedges, GeneList)

  }

  if(length(unique(all_pathways)) == 1) {
    return(unique(all_pathways))
  }

  n <- length(GeneList)
  new_matrix <- matrix(0, nrow = n, ncol = n)
  for (i in seq(1,n)) {
    for (j in seq(1,n)) {
      new_matrix[i, j] <- length(intersect(GeneList[[i]], GeneList[[j]]))/(length(c(GeneList[[i]], GeneList[[j]])))
    }
  }

  colnames(new_matrix) <- names(GeneList)
  rownames(new_matrix) <- names(GeneList)

  diag(new_matrix) <- 0
  new_matrix <- new_matrix/max(new_matrix)

  net <- igraph::graph_from_adjacency_matrix(new_matrix, weighted = TRUE, mode = 'undirected')

  net_leiden <- igraph::cluster_leiden(net, objective_function = 'modularity', resolution = resolution)

  cluster_groups <- split(net_leiden$names, net_leiden$membership)

  names_clusters <- lapply(cluster_groups, function(i) {

    all_words <- unlist(strsplit(i, ' ', fixed = TRUE))
    rm_idxs <- c(grep('go:', tolower(all_words)), grep('gse', tolower(all_words)), grep('kegg', tolower(all_words)))
    if(length(rm_idxs) > 0) {
      all_words <- all_words[-rm_idxs]
    }
    all_words <- gsub(',', '', all_words, fixed = TRUE)
    all_words <- gsub(';', '', all_words, fixed = TRUE)
    all_words <- gsub('.', '', all_words, fixed = TRUE)
    all_words <- gsub(':', '', all_words, fixed = TRUE)
    all_words <- gsub('(', '', all_words, fixed = TRUE)
    all_words <- gsub(')', '', all_words, fixed = TRUE)
    all_words <- all_words[which(!tolower(all_words) %in% c('of', 'and', 'by', 'in', 'to', 'via', 'sample', '-', '.', ';', '_', 'the', 'with'))]

    all_words <- all_words[which(!is_number(all_words))]
    all_words <- tolower(all_words)
    all_words <- unlist(lapply(all_words, SemNetCleaner::singularize))
    all_words[which(all_words == 'diseases')] <- 'disease'

    t_all_words <- sort(table(all_words), decreasing = TRUE)

    sel_words <- names(t_all_words)[which(t_all_words > stats::quantile(t_all_words, 0.9))]

    new_all_words <- i
    for(w in sel_words) {
      new_all_words <- new_all_words[-grep(w, new_all_words)]
    }

    new_all_words <- unlist(strsplit(new_all_words, ' ', fixed = TRUE))
    rm_idxs <- c(grep('go:', tolower(new_all_words)), grep('gse', tolower(new_all_words)), grep('kegg', tolower(new_all_words)))
    if(length(rm_idxs) > 0) {
      new_all_words <- new_all_words[-rm_idxs]
    }
    new_all_words <- gsub(',', '', new_all_words, fixed = TRUE)
    new_all_words <- gsub(';', '', new_all_words, fixed = TRUE)
    new_all_words <- gsub('.', '', new_all_words, fixed = TRUE)
    new_all_words <- gsub(':', '', new_all_words, fixed = TRUE)
    new_all_words <- gsub('(', '', new_all_words, fixed = TRUE)
    new_all_words <- gsub(')', '', new_all_words, fixed = TRUE)
    new_all_words <- new_all_words[which(!tolower(new_all_words) %in% c('of', 'and', 'by', 'in', 'to', 'via', 'sample', '-', '.', ';', '_', 'the', 'with'))]

    new_all_words <- new_all_words[which(!is_number(new_all_words))]
    new_all_words <- tolower(new_all_words)
    new_all_words <- unlist(lapply(new_all_words, SemNetCleaner::singularize))
    new_all_words[which(new_all_words == 'diseases')] <- 'disease'
    new_all_words <- c(new_all_words, all_words[which(!all_words %in% sel_words)])

    t_new_all_words <- sort(table(new_all_words), decreasing = TRUE)
    sel_new_words <- names(t_new_all_words)[which(t_new_all_words > stats::quantile(t_new_all_words, 0.9))]

    final_words <- c(sel_words, sel_new_words)

    if(length(final_words) == 0) {
      final_words <- names(t_all_words)[which(t_all_words == max(t_all_words))]
    }

    if('positive' %in% final_words & 'negative' %in% final_words) {
      final_words <- final_words[which(!final_words %in% c('positive', 'negative'))]
      final_words <- c(final_words, 'regulation')
    }

    if(length(unique(final_words)) > max.words) {
      sub_all_words <- all_words[which(all_words %in% final_words)]
      t_sub_all_words <- sort(table(sub_all_words), decreasing = TRUE)
      final_words <- names(t_sub_all_words)[seq(1,max.words)]
    }

    if(length(unique(final_words)) == 1) {
      sub_all_words <- all_words[which(!all_words %in% final_words)]
      t_sub_all_words <- sort(table(sub_all_words), decreasing = TRUE)
      if(length(t_sub_all_words) != 0) {
        if(!all(t_sub_all_words == min(t_sub_all_words))) {
          final_words <- c(final_words, names(t_sub_all_words)[1])
        }
      }
    }

    out <- paste0(unique(final_words), collapse = '.')

    return(out)

  })

  names(cluster_groups) <- names_clusters

  unique_groups <- unique(names(cluster_groups))

  cluster_groups <- lapply(unique_groups, function(name) {
    sub_groups <- cluster_groups[which(names(cluster_groups) == name)]
    names(sub_groups) <- NULL
    return(unique(unlist(sub_groups)))
  })

  names(cluster_groups) <- unique_groups

  unique_genelist <- unique(names(GeneList))
  united_GeneList <- lapply(unique_genelist, function(name) {
    sub_groups <- GeneList[which(names(GeneList) == name)]
    names(sub_groups) <- NULL
    return(unique(unlist(sub_groups)))
  })

  names(united_GeneList) <- unique_genelist

  group_avg_exp <- do.call(rbind, lapply(cluster_groups, function(anno) {
    genes <- unique(unlist(united_GeneList[anno]))
    if(length(genes) == 1) {
      return(SingleCellExperiment::logcounts(pseudo)[genes,])
    } else {
      return(Matrix::colMeans(SingleCellExperiment::logcounts(pseudo)[genes,]))
    }
  }))

  variable_padj <- do.call(cbind, lapply(fgseas_variable, function(i) {
    v <- -log10(i$padj)*sign(i$NES)
    names(v) <- i$pathway
    v <- v[fgseas_variable[[1]]$pathway]
    return(v)
  }))

  rownames(variable_padj) <- fgseas_variable[[1]]$pathway

  variable_padj <- apply(variable_padj, 1, function(i) {max(i, na.rm = TRUE)})
  variable_padj <- variable_padj[names(united_GeneList)]
  variable_padj[which(is.na(variable_padj))] <- 0

  if(!is.null(between_class_gsea)) {
    invariable_padj <- -log10(between_class_gsea$padj) * sign(between_class_gsea$NES)
    names(invariable_padj) <- between_class_gsea$pathway
    invariable_padj <- invariable_padj[names(united_GeneList)]
    invariable_padj[which(is.na(invariable_padj))] <- 0
  }

  groups_to_anno <- S4Vectors::DataFrame('Annotation' = unlist(cluster_groups),
                              'AnnotationGroup' = unlist(lapply(names(cluster_groups),
                                                                function(g) {rep(g,
                                                                                 length(cluster_groups[[g]]))})))

  groups_to_anno <- groups_to_anno[match(names(united_GeneList), groups_to_anno$Annotation),]

  avg_exp <- get_column_group_average(group_avg_exp, net_df$Stages)

  highest_exp_stage <- colnames(avg_exp)[apply(avg_exp, 1, which.max)]

  names(highest_exp_stage) <- rownames(avg_exp)

  significance_var <- as.vector(rowsum(abs(as.numeric(variable_padj)),
                                       groups_to_anno$AnnotationGroup)/as.matrix(table(groups_to_anno$AnnotationGroup)))

  if(!is.null(between_class_gsea)) {
    significance_invar <- as.vector(rowsum(abs(as.numeric(invariable_padj)),
                                           groups_to_anno$AnnotationGroup)/as.matrix(table(groups_to_anno$AnnotationGroup)))

    significance <- Matrix::rowMeans(cbind(significance_var,
                                   significance_invar))
  } else {
    significance <- significance_var
  }

  cluster_groups <- cluster_groups[unique(groups_to_anno$AnnotationGroup)]
  n <- length(cluster_groups)
  genes_matrix <- matrix(0, nrow = n, ncol = n)

  for (i in seq(1,n)) {
    for (j in seq(1,n)) {
      g1 <- unique(unlist(geneset_list[cluster_groups[[i]]]))
      g2 <- unique(unlist(geneset_list[cluster_groups[[j]]]))
      genes_matrix[i, j] <- length(intersect(g1, g2))/(length(c(g1, g2)))
    }
  }

  colnames(genes_matrix) <- unique(names(cluster_groups))
  rownames(genes_matrix) <- unique(names(cluster_groups))

  diag(genes_matrix) <- 0

  if(max(genes_matrix) != 0) {
    genes_matrix <- genes_matrix/max(genes_matrix, na.rm = TRUE)
  }

  network_geneset <- igraph::graph_from_adjacency_matrix(genes_matrix, weighted = TRUE, mode = 'undirected')

  igraph::V(network_geneset)$size <- (sqrt(significance))/max(sqrt(significance))*5
  igraph::V(network_geneset)$color <- stages_palette[highest_exp_stage[unique(groups_to_anno$AnnotationGroup)]]

  groups_to_anno$GroupStage <- highest_exp_stage[groups_to_anno$AnnotationGroup]

  groups_to_anno$pAdjTime <- variable_padj[groups_to_anno$Annotation]

  groups_to_anno$leadingEdges <- unlist(lapply(united_GeneList, function(i) {paste0(i, collapse = '-')}))

  top_genes_by_group <- lapply(unique(groups_to_anno$AnnotationGroup), function(i) {
    all_anno <- groups_to_anno$Annotation[which(groups_to_anno$AnnotationGroup == i)]
    all_genes <- unlist(united_GeneList[all_anno])

    top_genes <- sort(table(all_genes), decreasing = TRUE)

    out <- paste0(names(top_genes)[seq(1,min(10, length(top_genes)))], collapse = '-')

    return(out)

  })

  names(top_genes_by_group) <- unique(groups_to_anno$AnnotationGroup)

  top_genes_by_group <- unlist(top_genes_by_group)

  groups_to_anno$TopGenesByGroup <- top_genes_by_group[groups_to_anno$AnnotationGroup]

  rownames(groups_to_anno) <- groups_to_anno$Annotation

  plot_geneset_network(network_geneset)

  out <- S4Vectors::List('GeneSetDF' = groups_to_anno,
                         'Net' = network_geneset,
                         'GSEA' = fgseas_variable)

  if(!is.null(between_class_gsea)) {
    out$classGSEA <- between_class_gsea
  }

  return(out)

}
