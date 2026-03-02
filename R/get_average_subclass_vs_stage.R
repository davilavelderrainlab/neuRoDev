#' Get average expression divided into subclass and stages
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

  genes <- intersect(genes, rownames(net))
  if (!length(genes)) stop("No genes found in net rownames().")

  stages   <- net@colData[, "Stages"]
  subclass <- net@colData[, "SubClass"]

  stage_levels <- sort(unique(stages))
  sub_levels   <- sort(unique(subclass))

  out <- matrix(0, nrow = length(stage_levels), ncol = length(sub_levels),
                dimnames = list(stage_levels, sub_levels))

  X <- SummarizedExperiment::assays(net)[["logcounts"]]

  for (st in stage_levels) {
    idx <- which(stages == st)
    if (!length(idx)) next

    if (length(genes) == 1) {
      x <- as.numeric(X[genes, idx, drop = TRUE])
      s <- stats::sd(x)
      if (is.finite(s) && s > 0) {
        z <- (x - mean(x)) / s
      } else {
        z <- rep(0, length(x))
      }
      med <- tapply(z, subclass[idx], stats::median)
    } else {
      M <- as.matrix(X[genes, idx, drop = FALSE])
      mu <- matrixStats::rowMeans2(M)
      sdv <- matrixStats::rowSds(M)
      sdv[!is.finite(sdv) | sdv == 0] <- NA_real_
      M <- (M - mu) / sdv
      score <- colMeans(M, na.rm = TRUE)
      med <- tapply(score, subclass[idx], stats::median)
    }

    out[st, names(med)] <- med
    out[st, is.na(out[st, ])] <- 0
  }
  out
}
