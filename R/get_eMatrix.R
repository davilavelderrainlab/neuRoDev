#' Get the expression enrichment matrix
#'
#' @param net The reference network
#' @param genes Genes of which to compute the expression enrichment
#' @param nRand The number of random sampling to use (defaults to 100)
#'
#' @return A list with z-scores and p-values matrices divided into
#' subclass-stage pairs
#' @export
#'
#' @examples
#' m <- matrix(sample(seq(1,10, length.out=10000), 15000*100, replace = TRUE), ncol = 100)
#' rownames(m) <- paste0('Gene-', seq(1,15000))
#' colnames(m) <- paste0('Col-', seq(1,100))
#' net <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = m))
#' net$SubClass <- rep(c('A', 'B', 'C', 'D'), each = 25)
#' net$Stages <- rep(c('S1', 'S2', 'S3', 'S4'), each = 25)
#' emat <- get_eMatrix(net = net, genes = rownames(net)[seq(1,5)])
get_eMatrix <- function(net, genes, nRand=100) {

  Test <- get_average_subclass_vs_stage(net, genes)

  RandSets <- lapply(seq(1,nRand), function(i) {sample(rownames(net), length(genes))})

  RandSets <- lapply(RandSets, function(i) {get_average_subclass_vs_stage(genes=i, net = net)})

  Ave <- do.call("cbind",lapply(seq(1,ncol(RandSets[[1]])), function(j) {
    Matrix::rowMeans(do.call(cbind, lapply(RandSets, function(i) i[,j])))
  }))

  SDs <- do.call("cbind",lapply(seq(1,ncol(RandSets[[1]])), function(j) {
    apply(do.call(cbind, lapply(RandSets, function(i) i[,j])), 1, stats::sd)
  }))

  z <- (Test-Ave)/SDs
  p <- 2 * stats::pnorm(-abs(z))

  return(S4Vectors::List(z=z, p=p))

}
