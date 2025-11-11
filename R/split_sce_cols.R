#' Split a SingleCellExperiment object by column
#'
#' @param sce The SingleCellExperiment
#' @param colVec A vector that groups columns
#'
#' @return A list of SingleCellExperiment objects with one element per group
#'
#' @examples
#' m <- matrix(sample(seq(1,10, length.out=10000), 15000*100, replace = TRUE), ncol = 100)
#' rownames(m) <- paste0('Gene-', seq(1,15000))
#' colnames(m) <- paste0('Col-', seq(1,100))
#' net <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = m))
#' net$SubClass <- rep(c('A', 'B', 'C', 'D'), each = 25)
#' neuRoDev:::split_sce_cols(net, net$SubClass)
split_sce_cols <- function(sce, colVec) {
  groups <- as.character(sort(unique(colVec)))
  temp.sce.split <- lapply(groups, function(i) sce[,colVec==i])
  names(temp.sce.split) <- groups
  return(temp.sce.split)
}
