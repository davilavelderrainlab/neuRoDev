#' Gets the average by-column of columns belonging to a certain group (profiles)
#'
#' This function computes the column average of each column belonging to a given
#' group defined by the membership vector (profiles).
#'
#' @param M A matrix
#' @param group A membership vector for the columns of M
#' @param na.rm If NAs must be removed in the calculation. Defaults to TRUE
#'
#' @return A table with a single entry for each group with the profile.
#' @export
#'
#' @examples
#' get_column_group_average(M = matrix(seq(1,100), ncol=10),
#' group = c(rep(c('A','B','C'), each = 3), 'D'))
#' get_column_group_average(M = matrix(seq(1,100), ncol=10),
#' group = c(rep(c('A','B','C'), each = 3), 'D'))
get_column_group_average <- function(M,
                                     group,
                                     na.rm = TRUE) {

  group <- as.vector(group)

  M_sub <- M[,which(group %in% names(table(group))[which(table(group)>1)])]

  if(is.vector(M_sub)) {
    M_sub <- t(as.matrix(M_sub))
    return(t((rowsum(t(M_sub), colnames(M_sub)))/as.matrix(table(colnames(M_sub))))[1,])
  }

  M_sub_comp <- as.matrix(M[,which(!group %in% names(table(group))[which(table(group)>1)])])
  if(dim(M_sub_comp)[2] == 1) {
    colnames(M_sub_comp) <- names(table(group))[which(table(group)==1)]
  } else {
    colnames(M_sub_comp) <- group[match(colnames(M_sub_comp), colnames(M))]
  }

  group_sub <- group[which(group %in% names(table(group))[which(table(group)>1)])]

  out_sub <- do.call(cbind, lapply(split(seq_along(group_sub),
                                     group_sub),
                function(cols) {
                  if(length(cols) == 1) {
                    return(M_sub[,cols])
                  } else {
                    if(is.vector(M_sub)) {
                      return(mean(M_sub[cols], na.rm = na.rm))
                    } else {
                      return(Matrix::rowMeans(M_sub[,cols], na.rm = na.rm))
                    }
                  }
                 }))

  out <- cbind(M_sub_comp, out_sub)

  out <- out[,unique(group)]

  return(out)

}
