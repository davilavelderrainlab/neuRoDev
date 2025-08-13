#' Sum of each element in a group
#'
#' @param v A named numerical vector. The names define the groups
#' @param perc A boolean variable. If TRUE, the result is given as a percentage
#' out of the total numbers (sum between groups = 1)
#' @param sort A boolean variable. If TRUE, results will be sorted
#' (decreasing=TRUE)
#' @param reciprocal A boolean variable. If TRUE, the vector will be first turned
#' into 1/v
#' @param compute_means A boolean variable. If TRUE, it computes the mean instead
#' of the sum of values belonging to a given group.
#'
#' @return The sum of each element in a group, where the groups are the names of
#' the vector
#'
#' @examples
#' v <- seq(1,100)
#' names(v) <- rep(c('A', 'B', 'C', 'D'), each = 25)
#' neuRoDev:::value_table(v)
#' neuRoDev:::value_table(v, perc=TRUE)
#' neuRoDev:::value_table(v, perc=TRUE, reciprocal=TRUE)
value_table <- function(v,
                        perc=FALSE,
                        sort=TRUE,
                        reciprocal=FALSE,
                        compute_means=FALSE) {

  if(reciprocal) {

    tab <- c()
    for(i in unique(names(v))) {
      f_v <- v[which(names(v) == i)]
      if(compute_means) {
        tab <- c(tab, mean(1/f_v))
      } else {
        tab <- c(tab, sum(1/f_v))
      }
    }
    names(tab) <- unique(names(v))

  } else {
    tab <- c()
    for(i in unique(names(v))) {
      f_v <- v[which(names(v) == i)]
      if(compute_means) {
        tab <- c(tab, mean(f_v))
      } else {
        tab <- c(tab, sum(f_v))
      }
    }
    names(tab) <- unique(names(v))
  }

  if(perc) {
    tab <- tab / sum(tab)
  }

  if(sort) {
    tab <- sort(tab,
                decreasing = TRUE)
  }

  return(tab)
}
