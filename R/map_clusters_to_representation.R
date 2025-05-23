#' Wrap-up function to compute a series of summary information given an expression
#' matrix and a membership vector
#'
#' @param M An expression matrix
#' @param group A membership vector
#' @param BOGs A boolean variable to define if the BOGs will be computed
#'
#' @return A list with the profiles (get_column_group_average),
#' signatures (FC_signatures) and bag-of-genes
#' (run_wilcox_differential + filter_dge_out), and the differential expression
#' table (run_wilcox_differential)
#' @export
#'
#' @examples
#' map_clusters_to_representation(M = matrix(seq(1,100), ncol=10),
#' group = c(rep(c('A','B','C'), each = 3), 'D'))
map_clusters_to_representation <- function(M,
                                           group,
                                           BOGs=FALSE) {

  group_counts <- table(group)
  groups <- sort(unique(group))

  P <- get_column_group_average(M,
                                group)

  S <- FC_signatures(P)

  if(BOGs) {

    wilcox_dif_out <- lapply(groups, function(i) {
      run_wilcox_differential(M, binGroup = ifelse(group==i, 1, 0))
    })

    names(wilcox_dif_out) <- groups

    BOG <- lapply(wilcox_dif_out, function(i) {
      rownames(filter_dge_out(i))
    })

  } else {

    wilcox_dif_out <- NULL
    BOG <- NULL

  }

  out <- S4Vectors::List(group_counts=group_counts,
                         P=P,
                         S=S,
                         BOG=BOG,
                         DE_out=wilcox_dif_out)
  return(out)
}
