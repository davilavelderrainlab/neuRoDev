#' Get a matrix from a list
#'
#' @param l A list of named numerical vectors
#'
#' @return A matrix with as columns each name of any individual named vector and as
# rows all elements of the list
#'
#' @examples
#' v1 <- seq(2,12)
#' names(v1) <- c(rep('A', 5), rep('B', 6))
#' v2 <- seq(0,10)
#' names(v2) <- c(rep('B', 5), rep('C', 6))
#' l <- list('X' = v1, 'Y' = v2)
#' neuRoDev:::get_matrix_from_list(l)
get_matrix_from_list <- function(l) {

  all_names <- unlist(lapply(l, names))

  m <- matrix(0, ncol = length(unique(all_names)), nrow = length(l))
  colnames(m) <- unique(all_names)
  rownames(m) <- names(l)

  for(i in seq_len(length(l))) {
    m[i, names(value_table(l[[i]]))] <- value_table(l[[i]])
  }

  m <- m[,order(colnames(m))]

  return(m)
}
