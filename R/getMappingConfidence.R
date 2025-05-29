#' getMappingConfidence
#'
#' @inheritParams umap_plot_same_layout
#' @inheritParams reference_signatures_correlation
#' @param mapped_obj The full mapped object derived from umap_plot_same_layout or
#' add_to_reference (containing both the New and Original elements)
#'
#' @return The confidence and mappability of the given mapping
#' @export
#'
#' @examples
#' set.seed(123)
#' S <- FC_signatures(matrix(runif(200,0,10), ncol = 10))
#' rownames(S) <- paste0('Gene-', seq(1, dim(S)[1]))
#' refS <- FC_signatures(matrix(runif(200,0.1,7), ncol = 10))
#' colnames(refS) <- paste0('Reference-', seq(1, dim(refS)[2]))
#' rownames(refS) <- paste0('Gene-', seq(1, dim(refS)[1]))
#' annotated_M <- reference_signatures_correlation(S, refS)
#' new_profiles <- matrix(runif(80,0,10), ncol = 4)
#' rownames(new_profiles) <- paste0('Gene-', seq(1, dim(new_profiles)[1]))
#' colnames(new_profiles) <- paste0('New-', seq(1, dim(new_profiles)[2]))
#' new_clusterS <- FC_signatures(new_profiles)
#' rownames(new_clusterS) <- paste0('Gene-', seq(1, dim(new_clusterS)[1]))
#' colnames(new_clusterS) <- paste0('New-', seq(1, dim(new_clusterS)[2]))
#' new_M <- reference_signatures_correlation(new_clusterS, refS)
#' mapping <- add_to_reference(annotated_M,
#' new_M,
#' annotated_M$`Best.Assignment`)
#' getMappingConfidence(mapped_obj = mapping, signatures_cor = new_M)
getMappingConfidence <- function(mapped_obj,
                                 signatures_cor) {

  new_clusters <- rownames(signatures_cor)
  umap_obj <- mapped_obj$New$umap_obj
  old_umap_obj <- mapped_obj$Original$umap_obj

  sub_signatures_cor <- get_correlation_values(signatures_cor)

  sub_sub_signatures_cor <- sub_signatures_cor[new_clusters,]

  new_d_matrix <- igraph::distances(mapped_obj$New$umap_obj$refined_network,
                                    weights = (1-igraph::E(mapped_obj$New$umap_obj$refined_network)$weight))

  old_d_matrix <- igraph::distances(mapped_obj$Original$umap_obj$refined_network,
                                      weights = (1-igraph::E(mapped_obj$Original$umap_obj$refined_network)$weight))

  sub_new_d_matrix <- new_d_matrix[new_clusters, new_clusters]

  accuracy <- signatures_cor[new_clusters, "Best.Value"]

  ref_precision <- mean(Matrix::colMeans(old_d_matrix))

  min_precision <- min(igraph::E(mapped_obj$Original$umap_obj$refined_network)$weight)

  precision <- unlist(lapply(Matrix::colMeans(sub_new_d_matrix, na.rm = TRUE),
                             function(i) {
                                min(i, ref_precision)
                             }))

  precision <- min_max_normalize(c(ref_precision, precision, min_precision))

  precision <- precision[seq(2, (length(precision)-1))]

  score <- sqrt(accuracy*precision)

  return(S4Vectors::List(Confidence = S4Vectors::List(Accuracy = accuracy,
                                                      Precision = precision,
                                                      Score = score),
                         Mappability = S4Vectors::List(Accuracy = mean(accuracy),
                                                       Precision = mean(precision),
                                                       Score = mean(score)),
                         UMAP = mapped_obj$New$umap_plot))
}
