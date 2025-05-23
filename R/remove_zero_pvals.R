#' Remove pvalues equal to 0
#'
#' This function removes the pvalues equal to zeros and replaces them with a pvalue
#' equal to the minimum pvalue (excluding those equal to 0) times 0.1
#'
#' @param xVec A vector containing pvalues
#'
#' @return The vector with the modified 0 pvalues
#' @export
#'
#' @examples
#' remove_zero_pvals(c(rep(0, 5), seq(0.01,1,0.1)))
remove_zero_pvals <- function(xVec) {

  xVec[xVec==0] <- min(xVec[xVec!=0])*0.1
  return(xVec)
}
