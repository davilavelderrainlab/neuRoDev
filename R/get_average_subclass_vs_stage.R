#' Get average expression divided into subclass and signatures
#'
#' @param net The reference network
#' @param genes Genes of which to compute the average expression
#'
#' @return A matrix with the average expression by subclass and stage
#'
#' @examples
#' m <- matrix(sample(seq(1,10, length.out=10000), 15000*100, replace = TRUE), ncol = 100)
#' rownames(m) <- paste0('Gene-', seq(1,15000))
#' colnames(m) <- paste0('Col-', seq(1,100))
#' net <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = m))
#' net$SubClass <- rep(c('A', 'B', 'C', 'D'), each = 25)
#' net$Stages <- rep(c('S1', 'S2', 'S3', 'S4'), each = 25)
#' neuRoDev:::get_average_subclass_vs_stage(net = net, genes = rownames(m)[seq(1,3)])
get_average_subclass_vs_stage <- function(net, genes) {

  genes <- genes[which(genes %in% rownames(net))]

  networkByStage <- split_sce_cols(net, net@colData[,"Stages"])

  o <- lapply(networkByStage, function(i) {
    S <- SummarizedExperiment::assays(i)[['logcounts']][genes,]
    S <- t(scale(t(S)))
    return(split(Matrix::colMeans(S),i@colData[,"SubClass"]))
  }
  )
  o <- lapply(o, function(j) sapply(j, stats::median))
  groupVarOrder <- sort(unique(net@colData[,"SubClass"]))
  o <- do.call(rbind, lapply(o, function(i) i[groupVarOrder]))
  colnames(o) <- groupVarOrder
  o[is.na(o)] <- 0
  return(o)

}
