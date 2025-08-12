#' Network smoothing
#'
#' @param network The igraph network
#' @param scores The nodes scores to smooth
#' @param n_nearest_smooth The number of nearest neighbors to consider for
#' smoothing. Defaults to the number of nodes.
#'
#' @return The smoothed scores
#'
#' @examples
#' set.seed(123)
#' g <- igraph::erdos.renyi.game(n = 10,
#' p.or.m = 0.3,
#' type = "gnp",
#' directed = FALSE)
#' igraph::E(g)$weight <- sample(seq(1,10),
#' length(igraph::E(g)),
#' replace = TRUE)
#' scores <- stats::runif(igraph::vcount(g), min = 0, max = 100)
#' network_smoothing(g, scores)
network_smoothing <- function(network,
                              scores,
                              n_nearest_smooth=NULL) {

  if(is.null(n_nearest_smooth)) {
    n_nearest_smooth <- length(igraph::V(network)$name)
  }

  adj_matrix <- igraph::as_adjacency_matrix(network, attr = 'weight')

  diag(adj_matrix) <- 1

  similarities <- t(apply(adj_matrix, 1, function(r) {
    r[order(r, decreasing = TRUE)]
  }))

  indexes <- t(apply(adj_matrix, 1, function(r) {
    order(r, decreasing = TRUE)
  }))

  n_nearest_smooth <- min(n_nearest_smooth, dim(indexes)[2])

  smoothed_scores <- unlist(lapply(seq(1, nrow(similarities)), function(i) {

    sub_s <- similarities[i,]
    sub_i <- indexes[i,]

    closest_clusters <- rownames(similarities)[sub_i]

    sub_scores <- scores[sub_i]

    new_sub_s <- rep(0, length(sub_s))
    names(new_sub_s) <- names(sub_scores)
    new_sub_s[seq(1, n_nearest_smooth)] <- sub_s[seq(1, n_nearest_smooth)]

    new_weighted_scores <- stats::weighted.mean(x = sub_scores, w = new_sub_s)

  }))

  names(smoothed_scores) <- rownames(similarities)
  smoothed_scores <- smoothed_scores[names(scores)]


  return(smoothed_scores)


}
