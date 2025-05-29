#' Defines if a character is a number
#'
#' @param x A character
#'
#' @return TRUE if it's a number, FALSE otherwise
#' @export
#'
#' @examples
#' is_number('1')
is_number <- function(x) {
  grepl("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$", x)
}
