#' It returns a membership matrix for each element in the input list
#'
#' @param InList An input list
#'
#' @return A sparse membership matrix
#' @export
#'
#' @examples
#' list_to_assignment_matrix(list('A' = c(1,2,3),
#' 'B' = c(4,3,5),
#' 'C' = c(2,3,5)))
list_to_assignment_matrix <- function(InList) {

  Fs <- sort(unique(unlist(InList)))
  x <- do.call(cbind, lapply(InList, function(i) {
    as.numeric(Fs%in%i)
  }))

  rownames(x) <- Fs
  colnames(x) <- names(InList)
  return(methods::as(x, "sparseMatrix"))
}
