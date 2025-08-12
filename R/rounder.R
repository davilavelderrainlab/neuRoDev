#' Rounder
#'
#' @param n A number
#'
#' @return The higher round number
#'
#' @examples
#' rounder(1224.21)
rounder <- function(n) {

  n_split <-  unlist(strsplit(as.character(round(n)), ''))

  rounded <- as.numeric(paste0(c('1', rep(0, length(n_split))), collapse = ''))

  return(rounded)

}
