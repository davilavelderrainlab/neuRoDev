#' Get pseudobulk from an expression matrix
#'
#' @param exp_matrix The expression matrix (raw counts)
#' @param membership_vector The membership vector on which to group the
#' `exp_matrix` (e.g., clusters)
#'
#' @return A SingleCellExperiment object with the pseudobulk counts and normalized
#' logcounts
#' @export
#'
#' @examples
#' exp_matrix <- matrix(sample(seq(1,10, length.out=1000), 100000, replace = TRUE), ncol = 100)
#' membership_vector <- c(rep('A', 25), rep('B', 25), rep('C', 25), rep('D', 25))
#' pseudobulk <- get_pseudobulk(exp_matrix, membership_vector)
get_pseudobulk <- function(exp_matrix,
                           membership_vector) {

  unique_mv <- sort(names(which(table(membership_vector)>1)))

  pseudo_counts <- t(rowsum(t(exp_matrix), membership_vector))

  rownames(pseudo_counts) <- rownames(exp_matrix)

  norm_counts <- log2((t(t(pseudo_counts)/Matrix::colSums(pseudo_counts, na.rm = TRUE))*1000000) + 1)

  out <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=pseudo_counts,
                                                                logcounts=norm_counts))

  if(is.null(colnames(out))) {
    colnames(out) <- colnames(pseudo_counts)
  }

  out@metadata$cell_count <- table(membership_vector)

  return(out)

}
