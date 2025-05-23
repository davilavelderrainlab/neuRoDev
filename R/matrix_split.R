#' Splits a matrix in submatrices based on a membership vector
#'
#' @param M A matrix
#' @param group A membership vector for the columns of the matrix
#'
#' @return A list of submatrices
#' @export
#'
#' @examples
#' matrix_split(matrix(sample(seq(0,1),100, replace=TRUE), ncol=10),
#' group=c(rep(c(1,2,3), each = 3), 4))
matrix_split <- function(M,
                         group) {

  m_list <- list()
  for(i in unique(group)) {
    m_list[[i]] <- M[,which(group == i)]
  }

  return(m_list)

}
