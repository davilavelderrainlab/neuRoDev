#' Computes the Wilcoxon test on a group of an expression matrix against the others
#'
#' @param expMat An expression matrix
#' @param group A membership vector for the columns of expMat
#' @param qG A query group
#' @param orderOut A boolean variable to determine if the output is ordered
#' based on the AUC (decreasing=TRUE)
#'
#' @return The results of the Wilcoxon test: a dataframe of the differential
#' expression analysis with the main
#' statistics, like the p-value of each gene, the Bonferroni adjusted pvalue (padj),
#' the log Fold-Change, the AUC, the average expression (avgExpr), the percentage of
#' non zero values (pct_in), and a score which is the -log10 of the adjusted
#' pvalue multiplied by the sign of the logFC.
#' @export
#'
#' @examples
#' set.seed(123)
#' expMat = matrix(runif(200,0,10), ncol = 10)
#' rownames(expMat) <- paste0('Gene-', seq(1, nrow(expMat)))
#' get_pairwise_wilcoxauc(expMat = expMat,
#' group = c(rep(c(1,2,3), each = 3), 4), qG = 2)
get_pairwise_wilcoxauc <- function(expMat,
                                   group,
                                   qG,
                                   orderOut=TRUE) {

  if(!is.numeric(group)) {
    group <- make.names(group)
  }

  tGs <- unique(group)[!unique(group)%in%qG]

  if(is.null(rownames(expMat))) {
    rownames(expMat) <- paste0('Feature', seq(1, nrow(expMat)))
  }

  O <- lapply(tGs, function(i) {

    temp <- scran::pairwiseWilcox(cbind(expMat[,group==qG],
                                        expMat[,group==i]),
                                  ifelse(c(group[group==qG],
                                           group[group==i])==qG, 1, 0))
    temp <- temp[[1]][[1]]
    temp$padj <- stats::p.adjust(remove_zero_pvals(temp$p.value), "bonferroni")

    if(sum(group == i) > 1) {
      mean_i <- Matrix::rowMeans(expMat[,group==i])
      pct_out <- Matrix::rowMeans(expMat[,group==i] != 0)*100
    } else {
      mean_i <- expMat[,group==i]
      pct_out <- as.numeric(expMat[,group==i] != 0)*100
    }

    if(sum(group == qG) > 1) {
      mean_qG <- Matrix::rowMeans(expMat[,group==qG])
      pct_in <- Matrix::rowMeans(expMat[,group==qG] != 0)*100
    } else {
      mean_qG <- expMat[,group==qG]
      pct_in <- as.numeric(expMat[,group==qG] != 0)*100
    }

    logFC <- mean_qG-mean_i

    temp$logFC <- logFC

    temp$avgExpr <- mean_qG

    temp$score <- -log10(temp$padj) * ifelse(temp$logFC>0, 1, -1)

    temp$AUC <- 1-temp$AUC

    if(orderOut) {

      temp <- temp[order(temp$AUC, decreasing = TRUE),]

    }

    temp$auc <- temp$AUC
    temp$AUC <- NULL

    temp$pct_in <- pct_in
    temp$pct_out <- pct_out

    temp$feature <- rownames(temp)

    return(temp)

  })

  names(O) <- unique(tGs)

  return(O)
}
