#' Computes the UMAP given a matrix
#'
#' Given a matrix M, it computes the UMAP and returns the
#' object, with the KNN graph, clustered with either Louvain, Leiden or
#' Walktrap (defined with method)
#' @param M A matrix
#' @param method A method for clustering an igraph network
#' @param resolution A resolution parameter in case the method is leiden.
#' Defaults to 1
#' @param n_neighbors The number of neighbors to consider in the UMAP
#' @param transpose A boolean variable that defines if the M matrix has to
#' be trasposed. Default is TRUE, as the matrix returned by
#' reference_signatures_correlation has to be transposed.
#' @param n_components The number of components, or dimensions, to compute
#' @param refine_network A boolean variable (defaults to TRUE) that defines if
#' a correlation-distance based network has to be computed and visualized in the
#' 2D plot
#' @param weight_quantile The quantile below which the edges will be put to NA
#' @param weights_normalization_coef A normalization coefficient to reduce the
#' width of the edges in the plots. It defaults to the number of nodes of the
#' network (with a maximum of 1000)
#'
#' @return A list containing the UMAP object, the KNN network and the output
#' of the clustering
#' @export
#'
#' @examples
#' set.seed(123)
#' S <- FC_signatures(matrix(runif(200,0,10), ncol = 10))
#' rownames(S) <- paste0('Gene-', seq(1, dim(S)[1]))
#' refS <- FC_signatures(matrix(runif(200,0.1,7), ncol = 10))
#' colnames(refS) <- paste0('Reference-', seq(1, dim(refS)[2]))
#' rownames(refS) <- paste0('Gene-', seq(1, dim(refS)[1]))
#' M <- reference_signatures_correlation(S, refS)
#' umap_graph_clustering(M)
umap_graph_clustering <- function(M,
                                  method='louvain',
                                  resolution=NULL,
                                  n_neighbors=15,
                                  transpose=TRUE,
                                  n_components=2,
                                  refine_network=TRUE,
                                  weight_quantile=NULL,
                                  weights_normalization_coef=1) {

  M <- get_correlation_values(M)

  if(transpose) {M <- t(M)}

  if(dim(M)[1]>dim(M)[2]) {

    n_neighbors <- min(n_neighbors, dim(M)[2]-1)

    Coords <- umap::umap(t(M),
                         n_neighbors=n_neighbors,
                         n_components=n_components)
    each_n <- n_neighbors
  } else {
    n_neighbors <- min(n_neighbors, dim(M)[1])
    Coords <- umap::umap(t(M),
                         n_neighbors=n_neighbors,
                         n_components=n_components)
    each_n <- n_neighbors
  }

  indexes <- Coords$knn$indexes
  distances <- Coords$knn$distances

  TE <- as.numeric(t(indexes))
  Ds <- as.numeric(t(distances))
  EL <- cbind(rep(rownames(indexes),
                  each=each_n),
              rownames(indexes)[TE])

  G <- igraph::graph_from_edgelist(EL,
                                   directed = FALSE)
  igraph::E(G)$weight <- Ds
  G <- igraph::simplify(G)

  network_to_use <- G

  if(refine_network) {
    snn_network <- compute_snn_network(signatures_cor = t(M),
                                       n_neighbors = n_neighbors)
    network <- snn_network$network
    edges <- build_edges_df(network,
                            layout=Coords$layout,
                            weight_quantile = weight_quantile,
                            weights_normalization_coef=weights_normalization_coef)

    network_to_use <- network

    distances_network <- 1-snn_network$`sNN.Matrix`

    distances <- t(apply(distances_network, 1, function(r) {
      r[order(r)]
    }))
    indexes <- t(apply(distances_network, 1, function(r) {
      order(r)
    }))

    Coords$refined_network$indexes <- indexes
    Coords$refined_network$distances <- distances
    Coords$refined_network$edges <- edges

  } else {
    network <- NULL
    Coords$refined_network <- network
  }

  if(tolower(method) == 'louvain') {
    CC <- igraph::cluster_louvain(network_to_use)
  }
  if(tolower(method) == 'leiden') {
    if(is.null(resolution)) {
      resolution <- 1
    }
    CC <- igraph::cluster_leiden(network_to_use,
                                 resolution_parameter = resolution)
  }
  if(tolower(method) == 'walktrap') {
    CC <- igraph::cluster_walktrap(network_to_use)
  }

  out <- S4Vectors::List(umap_out=Coords,
                         umap_knn_igraph=G,
                         refined_network=network,
                         clustering_out=CC,
                         n_components=n_components,
                         n_neighbors=n_neighbors)

  if(refine_network) {
   out[['snn_network']] <- snn_network
  }

  return(out)
}
