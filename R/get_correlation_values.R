#' Get correlation values from a correlation_df object
#'
#' @param correlation_df A correlation dataframe obtained from
#' `reference_signatures_correlation`
#'
#' @return Only the correlation values as a matrix
#' @export
#'
#' @examples
#' set.seed(123)
#' S <- FC_signatures(matrix(runif(200,0,10), ncol = 10))
#' rownames(S) <- paste0('Gene-', seq(1, dim(S)[1]))
#' refS <- FC_signatures(matrix(runif(200,0.1,7), ncol = 10))
#' colnames(refS) <- paste0('Reference-', seq(1, dim(refS)[2]))
#' rownames(refS) <- paste0('Gene-', seq(1, dim(refS)[1]))
#' cor_df <- reference_signatures_correlation(S, refS)
#' get_correlation_values(cor_df)
get_correlation_values <- function(correlation_df) {

  if(any(endsWith(tolower(colnames(correlation_df)), 'value'))) {
    sub_correlation_df <- correlation_df[,seq((max(which(endsWith(tolower(colnames(correlation_df)), 'value'))))+1,
                                              dim(correlation_df)[2])]
  } else {
    sub_correlation_df <- correlation_df
  }

  to_remove_idxs <- which(apply(sub_correlation_df, 2, function(i) {mean(is_number(i)) != 1}))

  if(length(to_remove_idxs) > 0) {
    to_remove_idxs <- seq(min(to_remove_idxs), dim(sub_correlation_df)[2])
    sub_correlation_df <- sub_correlation_df[,-to_remove_idxs]
  }

  sub_correlation_df <- as.matrix(sub_correlation_df)

  return(sub_correlation_df)

}
