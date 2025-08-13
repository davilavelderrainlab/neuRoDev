#' Min-Max normalization
#'
#' @param vector A numerical vector
#'
#' @return The normalized numerical vector with a new distribution between 0 and
#' 1
#'
#' @examples
#' v <- c(1,6,5,3,7,8,3,2)
#' n_v <- neuRoDev:::min_max_normalize(v)
min_max_normalize <- function(vector) {
  norm <- (vector - min(vector))/(max(vector) - min(vector))
  return(norm)
}
