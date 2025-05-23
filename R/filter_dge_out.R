#' Filters genes based on differential expression parameters
#'
#' @param degOut A differential expression table as defined by run_wilcox_differential
#' @param cutFrac A threshold for the percentage of non-zero fature values
#' @param cutFC A threshold for the Fold-Change
#' @param cutPadj A threshold for the adjusted pvalue
#' @param orderOut A boolean variable to order the result based on the AUC (decreasing=TRUE)
#'
#' @return The differential expression table filtered
#' @export
#'
#' @examples
#' filter_dge_out(run_wilcox_differential(expMat = matrix(seq(1,200), ncol=10),
#' binGroup = c(rep(c(1,2,3),
#' each = 3),
#' 4)))
filter_dge_out <- function(degOut,
                           cutFrac=0.25,
                           cutFC=0.05,
                           cutPadj=0.01,
                           orderOut=TRUE) {

  o <- degOut
  o <- o[o$pct_in>cutFrac,]
  o <- o[o$padj<cutPadj,]
  o <- o[o$logFC>cutFC,]

  if(orderOut) {
    o <- o[order(o$auc, decreasing = TRUE),]
  }

  return(o)
}
