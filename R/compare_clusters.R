#' Find Statistically Significant clusters
#'
#' @inheritParams plot_eTrace
#' @param return_tests A boolean variable to define if the results of the
#' Wilcoxon tests should be returned or not. Defaults to TRUE.
#'
#' @return Indexes of clusters that are statistically significant and the result
#' of the statistical tests (if `return_tests` is TRUE)
#' @export
#'
#' @examples
#' m <- matrix(sample(seq(1,10, length.out=10000), 15000*100, replace = TRUE), ncol = 100)
#' rownames(m) <- paste0('Gene-', seq(1,15000))
#' colnames(m) <- paste0('Col-', seq(1,100))
#' net <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = m))
#' net$SubClass <- rep(c('A', 'B', 'C', 'D'), each = 25)
#' subclass_palette <- c('A' = 'red', 'B' = 'blue', 'C' = 'green', 'D' = 'yellow')
#' net$SubClass_color <- subclass_palette[net$SubClass]
#' net$Stages <- rep(c('S1', 'S2', 'S3', 'S4'), each = 25)
#' stages_palette <- c('S1' = 'pink', 'S2' = 'orange', 'S3' = 'violet', 'S4' = 'black')
#' net$Stages_color <- stages_palette[net$Stages]
#' net$X_coord <- sample(seq(1,2, length.out = 1000), size = ncol(net), replace = TRUE)
#' net$Y_coord <- sample(seq(1,2, length.out = 1000), size = ncol(net), replace = TRUE)
#' edges_from <- sample(colnames(net), size = 200, replace = TRUE)
#' edges_to <- sample(colnames(net), size = 200, replace = TRUE)
#' edges_from_x <- net$X_coord[match(edges_from, colnames(net))]
#' edges_from_y <- net$Y_coord[match(edges_from, colnames(net))]
#' edges_to_x <- net$X_coord[match(edges_to, colnames(net))]
#' edges_to_y <- net$Y_coord[match(edges_to, colnames(net))]
#' edges_weight <- sample(seq(0,1, length.out=1000), length(edges_from), replace = TRUE)
#' edges_df <- data.frame('from' = edges_from,
#' 'to' = edges_to,
#' 'weight' = edges_weight,
#' 'from.x' = edges_from_x,
#' 'from.y' = edges_from_y,
#' 'to.x' = edges_to_x,
#' 'to.y' = edges_to_y)
#' net@metadata$network$edges <- edges_df
#' SummarizedExperiment::rowData(net)$informative <- sample(c(TRUE, FALSE), size = nrow(net),
#' replace = TRUE)
#' compare_clusters(net, genes = 'Gene-1')
compare_clusters <- function(net,
                             genes,
                             expression_enrichment = TRUE,
                             pval_threshold = 0.05,
                             return_tests = TRUE) {

  if(expression_enrichment) {
    v <- get_eTrace(net, genes)$z
  } else {
    v <- Matrix::colMeans(SingleCellExperiment::logcounts(net)[genes,,drop=FALSE])
  }
  fit <- stats::lm(v ~ 0+interaction(net$SubClass,net$Stages))
  summ <- summary(fit)
  summ <- as.data.frame(summ$coefficients)
  rownames(summ) <- gsub('interaction(net$SubClass, net$Stages)', '', rownames(summ), fixed = TRUE)
  which_sig_overall <- rownames(summ)[which(summ$`t value` > 0 & stats::p.adjust(summ$`Pr(>|t|)`, 'fdr') < pval_threshold)]
  sig_stages <- unlist(lapply(strsplit(which_sig_overall, '.', fixed = TRUE), function(i) {i[2]}))
  sig_subclasses <- unlist(lapply(strsplit(which_sig_overall, '.', fixed = TRUE), function(i) {i[1]}))

  by_stage <- lapply(unique(sig_stages), function(stage) {
    sub_net <- net[,which(net$Stages == stage)]
    if(length(unique(sub_net$SubClass))==1) {
      return(paste0(unique(sub_net$SubClass), '.', stage))
    }
    score <- v[which(net$Stages == stage)]
    p <- lapply(intersect(unique(sub_net$SubClass),sig_subclasses), function(subclass) {
      stats::wilcox.test(score[which(sub_net$SubClass == subclass)], score[which(sub_net$SubClass != subclass)], alternative = 'greater')
    })
    names(p) <- paste0(intersect(unique(sub_net$SubClass),sig_subclasses), '.', stage)
    return(p)
  })

  by_stage_pval <- stats::p.adjust(unlist(lapply(by_stage, function(i) {lapply(i, function(x) x$p.value)})), 'fdr')
  by_stage_pval <- intersect(names(by_stage_pval)[which(by_stage_pval < pval_threshold)], which_sig_overall)
  idxs <- which(interaction(net$SubClass, net$Stages)%in%by_stage_pval)

  if(!return_tests) {
    return(idxs)
  } else {
    names(by_stage) <- unique(sig_stages)
    return(S4Vectors::List('GlobalTest'=summ, 'ByStageTests'=by_stage, 'SignificantIdxs'=idxs))
  }

  return(idxs)

}

