#' Wilcoxon differential expression analysis
#'
#' This function computes the Wilcoxon test of genetic expression matrices given
#' a membership vector.
#'
#' @param expMat A genetic expression matrix
#' @param binGroup A numerical membership vector for the columns of expMat
#' @param orderOut A boolean variable to order the output based on the AUC
#' (decreasing = TRUE)
#'
#' @return It returns a dataframe of the differential expression analysis with
#' the main statistics, like the p-value of each gene, the Bonferroni adjusted
#' pvalue (padj), the log Fold-Change, the AUC, the average expression (avgExpr),
#' and a score which is the -log10 of the adjusted
#' pvalue multiplied by the sign of the logFC.
#' @export
#'
#' @examples
#' run_wilcox_differential(expMat = matrix(seq(1,200), ncol=10),
#' binGroup = c(rep(c(1,2,3), each = 3), 4))
run_wilcox_differential <- function(expMat,
                                    binGroup,
                                    orderOut=FALSE) {

  if(!is.numeric(binGroup)) {
    binGroup <- make.names(binGroup)
  }

  if(is.null(rownames(expMat))) {
    rownames(expMat) <- paste0('Feature', seq(1, nrow(expMat)))
  }
  binGroup[which(binGroup != binGroup[1])] <- 'Other'
  o <- scran::pairwiseWilcox(expMat, binGroup)

  o <- o[[1]][[1]]
  o$padj <- stats::p.adjust(remove_zero_pvals(o$p.value), "bonferroni")

  if(sum(binGroup == binGroup[1]) > 1) {
    mean_one <- Matrix::rowMeans(expMat[,binGroup==binGroup[1]])
    pct_in <- Matrix::rowMeans(expMat[,binGroup==binGroup[1]] != 0)*100
  } else {
    mean_one <- expMat[,binGroup==binGroup[1]]
    pct_in <- as.numeric(expMat[,binGroup==binGroup[1]] != 0)*100
  }

  if(sum(binGroup != binGroup[1]) > 1) {
    mean_other <- Matrix::rowMeans(expMat[,binGroup!=binGroup[1]])
    pct_out <- Matrix::rowMeans(expMat[,binGroup!=binGroup[1]] != 0)*100
  } else {
    mean_other <- expMat[,binGroup!=binGroup[1]]
    pct_out <- as.numeric(expMat[,binGroup!=binGroup[1]] != 0)*100
  }

  logFC <- mean_one-mean_other

  o$logFC <- logFC

  o$avgExpr <- mean_one

  o$score <- -log10(o$padj) * ifelse(o$logFC>0, 1, -1)

  if(orderOut) {

    o <- o[order(o$AUC, decreasing = TRUE),]

  }

  o$auc <- o$AUC
  o$AUC <- NULL

  o$pct_in <- pct_in
  o$pct_out <- pct_out

  o$feature <- rownames(o)

  return(o)
}
