#' UMAP plot of the signatures
#'
#'  Given a matrix of correlation with signatures_cor,
#'  a color_attr parameter, a label attribute label_attr for the plot, it returns
#'  a UMAP object and a UMAP plot. The UMAP object can be given to the
#'  function with umap_obj. The clustering method of the UMAP is Louvain as
#'  default, but it can be also Walktrap ('walktrap') and Leiden ('leiden').
#'  A title for the UMAP can be given. A col_vector palette can be defined, one
#'  color for each color_attr. A legend can be added (legend=TRUE). A resolution
#'  parameter can be defined, if method='leiden'. The number of neighbors for
#'  umap_graph_clustering is defined by n_neighbors.
#' @inheritParams build_edges_df
#' @param signatures_cor The signature correlation matrix as defined in
#' reference_signatures_correlation (but only the correlation values)
#' @param refine_network A boolean variable (defaults to TRUE) that defines if
#' a correlation-distance based network has to be computed and visualized in the
#' 2D plot
#' @param color_attr The groups to show in the UMAP. It can be a column in
#' signatures_cor or a vector. It should match the names of col_vector, if given.
#' @param label_attr The labels to be used in the UMAP. If NULL, the groups given
#' in color_attr are shown. If no labels are wanted, set no_label to TRUE.
#' @param col_vector A named color vector palette with a color for each color_attr
#' @param umap_obj An already computed umap object as returned by
#' umap_graph_clustering
#' @param title The title of the UMAP
#' @param method The method to use in umap_graph_clustering. Can be either
#' 'leiden', 'louvain' (default) or 'walktrap'
#' @param legend A boolean variable that indicates if a legend has to be added
#' @param legend_name The name of the legend
#' @param resolution Resolution parameter for method='leiden'. Defaults to 1
#' @param n_neighbors The number of neighbors to be used in umap_graph_clustering.
#' It defaults to `min(15 + dim(signatures_cor)[1], dim(signatures_cor)[1])`
#' @param n_components The number of dimensions to use in the plot
#' @param show_edges A boolean variable to determine if the edges of the network
#' behind the plot should be shown
#' @param alpha The alpha vector for the colors in the UMAP. Defaults to 1
#' @param no_label A boolean variable. If TRUE, no label is shown.
#'
#' @return A list with the UMAP object and the UMAP plot
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
#' color_attr = M$`Best.Assignment`
#' umap_signature_plot(M,
#' color_attr = color_attr,
#' label_attr = color_attr)
umap_signature_plot <- function(signatures_cor,
                                refine_network=TRUE,
                                color_attr=NULL,
                                label_attr=NULL,
                                col_vector=NULL,
                                umap_obj=NULL,
                                title=NULL,
                                method='louvain',
                                legend=FALSE,
                                legend_name=NULL,
                                resolution=NULL,
                                n_neighbors=15,
                                n_components=2,
                                weight_quantile=NULL,
                                weights_normalization_coef=1,
                                show_edges=TRUE,
                                alpha=NULL,
                                no_label=FALSE) {

  if('umap_obj' %in% names(umap_obj)) {
    umap_obj <- umap_obj$umap_obj
  }

  if(is.null(label_attr) & !no_label) {
    label_attr <- color_attr
  }

  if(!is.null(label_attr) & no_label) {
    label_attr <- NULL
  }

  if(is.character(label_attr) && length(label_attr) == 1 && label_attr %in% colnames(signatures_cor)) {
    label_attr <- signatures_cor[,label_attr]
  }

  if(is.character(color_attr) && length(color_attr) == 1 && color_attr %in% colnames(signatures_cor)) {
    color_attr <- signatures_cor[,color_attr]
  }

  if(!is.null(color_attr) & is.null(col_vector)) {
    col_vector <- Polychrome::createPalette(length(unique(color_attr))+2,
                                c("#ffffff", "#ff0000", "#00ff00", "#0000ff", "#000000"))
    col_vector <- col_vector[seq(3,length(col_vector))]
    col_vector <- col_vector[seq(1,length(unique(color_attr)))]
    names(col_vector) <- unique(color_attr)
  }

  if(is.null(umap_obj)) {
    umap_obj <- umap_graph_clustering(signatures_cor,
                                      method = method,
                                      refine_network = refine_network,
                                      n_neighbors = n_neighbors,
                                      n_components = n_components,
                                      weight_quantile = weight_quantile,
                                      weights_normalization_coef = weights_normalization_coef)
  }

  layout <- umap_obj$umap_out$layout

  if(dim(layout)[2] != n_components) {
    stop('The given layout has a different number of dimensions to those wanted by the user')
  }

  network <- umap_obj$refined_network

  if(refine_network) {
    if(is.null(network)) {
      network <- compute_correlation_network(signatures_cor)
      edges <- build_edges_df(network,
                              layout=layout,
                              weight_quantile=weight_quantile,
                              weights_normalization_coef=weights_normalization_coef)
    } else {
      edges <- umap_obj$umap_out$refined_network$edges
    }
  } else {
    network <- NULL
    edges <- NULL
  }

  if(is.null(color_attr)) {
    color_attr <- umap_obj$clustering_out$membership
    names(color_attr) <- umap_obj$clustering_out$names
    color_attr <- color_attr[rownames(umap_obj$umap_out$layout)]

    col_vector <- Polychrome::createPalette(length(unique(color_attr))+2,
                                c("#ff0000", "#ffffff", "#00ff00", "#0000ff", "#000000"))
    col_vector <- col_vector[seq(3,length(col_vector))]
    col_vector <- col_vector[seq(1,length(unique(color_attr)))]
    names(col_vector) <- unique(color_attr)
  }

  if(n_components == 2) {
    layout <- cbind(layout, 'Size' = max(2,100/dim(layout)[1]))
  } else {
    layout <- cbind(layout, 'Size' = max(200,10000/dim(layout)[1]))
  }

  layout <- cbind(layout, 'Colors' = col_vector[as.character(color_attr)])

  layout <- cbind(layout, 'Group' = color_attr)

  if(!is.null(label_attr)) {
    layout <- cbind(layout, 'SubGroup' = label_attr)
  }

  layout <- layout[order(match(layout[,'Group'], names(col_vector))),]

  if(n_components == 2) {
    if(!is.null(label_attr)) {
      p <- umap_pointsize(layout,
                          color_attr = layout[,'Group'],
                          label_attr = layout[,'SubGroup'],
                          title = title,
                          legend = legend,
                          legend_name = legend_name,
                          only_new_points = FALSE,
                          edges = edges,
                          show_edges = show_edges,
                          alpha = alpha)
    } else {
      p <- umap_pointsize(layout,
                          color_attr = layout[,'Group'],
                          label_attr = NULL,
                          title = title,
                          legend = legend,
                          legend_name = legend_name,
                          only_new_points = FALSE,
                          edges = edges,
                          show_edges = show_edges,
                          alpha = alpha)
    }
  } else {
    p <- plotly::plot_ly(x = layout[,1],
                         y = layout[,2],
                         z = layout[,3],
                         type="scatter3d",
                         mode="markers",
                         colors = col_vector,
                         color = layout[,'Group'],
                         size = I(as.numeric(layout[,'Size'])))

    p$x$layoutAttrs <- list('Params' = list('title' = title,
                                            'scene' = list('xaxis' = list('title' = 'UMAP-1'),
                                                           'yaxis' = list('title' = 'UMAP-2'),
                                                           'zaxis' = list('title' = 'UMAP-3'))))
  }

  out <- S4Vectors::List(umap_obj = umap_obj,
                         umap_plot = p)

  return(out)

}
