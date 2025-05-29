#' Plots the UMAP given new clusters, keeping the same layout before and after
#' adding them
#'
#' @inheritParams umap_signature_plot
#' @param reference_df The reference correlation matrix.
#' @param new_clusters The new clusters names.
#' @param new_name The name/label that has to be given to new clusters.
#' @param new_points_col The color of the new clusters in the UMAP plot
#' @param only_new_points A boolean variable that defines if only the new clusters
#' are labelled
#' @param n_components The number of dimensions to use in the plot
#' @param show_edges A boolean variable to define if the plot should show edges
#' of the network behind the plot
#' @param n_increase The increase in size of new edges
#' @param mapping_quality A boolean to define if the mapping quality should be
#' computed
#' @param annotate A boolean to define if the nearest
#' neighbor annotation has to be computed
#' @param ... Additional parameters for `nearest_neighbor_annotation`
#'
#' @return The umap object and plot before and after adding the new clusters,
#' plus mapping quality metrics and nearest neighbor annotation
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
#' new_clusterS <- FC_signatures(matrix(runif(80,0,10), ncol = 4))
#' rownames(new_clusterS) <- paste0('Gene-', seq(1, dim(new_clusterS)[1]))
#' colnames(new_clusterS) <- paste0('New-', seq(1, dim(new_clusterS)[2]))
#' new_M <- reference_signatures_correlation(new_clusterS, refS)
#' umap_plot_same_layout(reference_df = annotated_M,
#' signatures_cor = new_M,
#' color_attr = annotated_M$`Best.Assignment`)
umap_plot_same_layout <- function(reference_df,
                                  signatures_cor,
                                  color_attr,
                                  mapping_quality=TRUE,
                                  annotate=TRUE,
                                  refine_network=TRUE,
                                  umap_obj=NULL,
                                  method="louvain",
                                  new_clusters=NULL,
                                  col_vector=NULL,
                                  new_name=NULL,
                                  label_attr=NULL,
                                  title=NULL,
                                  new_points_col="#FF0000",
                                  legend=FALSE,
                                  legend_name=NULL,
                                  only_new_points=FALSE,
                                  resolution=NULL,
                                  n_neighbors=15,
                                  n_components=2,
                                  weight_quantile=NULL,
                                  weights_normalization_coef=1,
                                  show_edges=TRUE,
                                  n_increase=0.2,
                                  alpha=NULL,
                                  no_label=FALSE,
                                  ...) {

  if('umap_obj' %in% names(umap_obj)) {
    umap_obj <- umap_obj$umap_obj
  }

  sub_signatures_cor <- get_correlation_values(signatures_cor)

  sub_reference_df <- get_correlation_values(reference_df)

  sub_reference_df <- sub_reference_df[,colnames(sub_signatures_cor)]

  if(is.null(label_attr) & !no_label) {
    label_attr <- color_attr
  }

  if(no_label) {
    label_attr <- NULL
  }

  if(is.character(label_attr) && length(label_attr) == 1 && label_attr %in% colnames(reference_df)) {
    label_attr <- reference_df[,label_attr]
  }

  if(is.character(color_attr) && length(color_attr) == 1 && color_attr %in% colnames(reference_df)) {
    color_attr <- reference_df[,color_attr]
  }

  if(is.null(new_clusters)) {
    new_clusters <- rownames(signatures_cor)
  }

  if(is.null(new_name)) {
    new_name <- new_clusters
  }

  if(is.null(names(new_points_col))) {
    if(length(new_points_col) == 1) {
      new_points_col <- rep(new_points_col, length(new_clusters))
    }
    names(new_points_col) <- new_clusters
  }

  if(length(color_attr) == dim(reference_df)[1]) {
    old_color_attr <- color_attr
    if(length(new_name) == 1) {
      color_attr <- c(color_attr, rep(new_name, dim(signatures_cor)[1]))
    } else {
      color_attr <- c(color_attr, new_name)
    }
  } else {
    old_color_attr <- color_attr[which(!color_attr %in% new_name)]
  }

  if(!is.null(label_attr)) {
    if(length(label_attr) == dim(reference_df)[1]) {
      old_label_attr <- label_attr
      if(length(new_name) == 1) {
        label_attr <- c(label_attr, rep(new_name, dim(signatures_cor)[1]))
      } else {
        label_attr <- c(label_attr, new_name)
      }
    } else {
      old_label_attr <- label_attr[which(!label_attr %in% new_name)]
    }
  } else {
    old_label_attr <- NULL
  }

  if(any(new_clusters %in% c(old_label_attr, old_color_attr))) {
    new_clusters[which(new_clusters %in% c(old_label_attr, old_color_attr))] <- paste0('New-', new_clusters[which(new_clusters %in% c(old_label_attr, old_color_attr))])
  }

  rownames(signatures_cor) <- new_clusters

  if(is.null(reference_df) && is.null(umap_obj)) {
    stop('Either a layout or a starting matrix need to be provided')
  }

  if(is.null(umap_obj)) {
    umap_obj <- umap_graph_clustering(sub_reference_df,
                                     method = method,
                                     refine_network = refine_network,
                                     n_neighbors = n_neighbors,
                                     n_components = n_components)
  }

  layout <- umap_obj$umap_out$layout

  if(dim(layout)[2] != n_components) {
    stop('The given layout has a different number of dimensions to those wanted by the user.
            The fewest dimensions will be considered, but taking the first dimensions of a layout
            with an higher number of dimensions not always is the right choice. Consider using a different
            layout, with the wanted dimensions.')
  }

  if(n_components > 2) {
    layout <- layout[,seq(1,3)]
  } else {
    layout <- layout[,seq(1,2)]
  }

  network <- umap_obj$refined_network

  if(refine_network) {
    if(is.null(network)) {
      network <- compute_snn_network(sub_reference_df,
                                     n_neighbors=n_neighbors)$network
      edges <- build_edges_df(network,
                              layout = layout,
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
  }

  if(length(col_vector) == length(unique(old_color_attr))) {
    if (!is.null(names(col_vector))) {
      old_names <- names(col_vector)
      col_vector <- c(col_vector, new_points_col[new_clusters])
      names(col_vector) <- c(old_names, new_clusters)
    } else {
      col_vector <- c(col_vector, new_points_col[new_clusters])
      names(col_vector) <- c(unique(old_color_attr), new_clusters)
    }
  }

  if(is.null(col_vector)) {
    col_vector <- Polychrome::createPalette(length(unique(old_color_attr))+2,
                                c("#ffffff", "#ff0000", "#00ff00", "#0000ff", "#000000"))
    col_vector <- col_vector[seq(3,length(col_vector))]
    col_vector <- col_vector[seq(1,length(unique(old_color_attr)))]
    if(length(new_points_col) == 1) {
      col_vector <- c(col_vector, rep(new_points_col, length(new_name)))
    } else {
      col_vector <- c(col_vector, new_points_col)
    }
    names(col_vector) <- unique(color_attr)
    col_vector <- col_vector[which(!is.na(names(col_vector)))]
  }

  combined_matrix <- rbind(S4Vectors::DataFrame(sub_reference_df),
                           S4Vectors::DataFrame(sub_signatures_cor))

  if(is.null(n_neighbors)) {n_neighbors <- min(15 + dim(sub_signatures_cor)[1],
                                               dim(combined_matrix)[2])
  }

  new_umap_obj <- umap_graph_clustering(M = combined_matrix,
                                        method = method,
                                        refine_network = refine_network,
                                        n_neighbors = n_neighbors,
                                        n_components = n_components,
                                        weight_quantile = weight_quantile,
                                        weights_normalization_coef = weights_normalization_coef)

  indexes <- new_umap_obj$umap_out$knn$indexes
  distances <- new_umap_obj$umap_out$knn$distances

  new_network <- new_umap_obj$refined_network

  if(refine_network) {

    new_edges <- new_umap_obj$umap_out$refined_network$edges
    distances <- new_umap_obj$umap_out$refined_network$distances
    indexes <- new_umap_obj$umap_out$refined_network$indexes

  } else {
    new_edges <- NULL
  }

  neighbors_x <- c()
  neighbors_y <- c()
  if(n_components > 2) {
    neighbors_z <- c()
  }

  distances_raw <- stats::dist(combined_matrix)
  distances_raw <- Matrix::forceSymmetric(as.matrix(distances_raw), uplo="L")
  distances <- t(apply(distances_raw, 1, function(r) {
    r[order(r)]
  }))
  indexes <- t(apply(distances_raw, 1, function(r) {
    order(r)
  }))

  for(c in new_clusters) {

    f_ind <- indexes[c,seq(2,dim(indexes)[2])]
    f_distances <- distances[c,seq(2,dim(indexes)[2])]

    to_remove <- which(rownames(indexes)[f_ind] %in% new_clusters)

    if(dplyr::between(length(to_remove), 1, ((dim(indexes)[2]-1) * 50 / 100)+1)) {
      f_ind <- f_ind[-to_remove]
      f_distances <- f_distances[-to_remove]
    }

    if(length(to_remove) > ((dim(indexes)[2]-1) * 50 / 100)) {
      new_combined_matrix <- combined_matrix[which(!rownames(combined_matrix) %in% new_clusters),]
      new_combined_matrix <- rbind(new_combined_matrix, combined_matrix[c,])

      if(is.null(n_neighbors)) {n_neighbors <- min(15 + length(new_clusters),
                                                   dim(new_combined_matrix)[2])
      }

      partial_umap_obj <- umap_graph_clustering(new_combined_matrix,
                                                method = method,
                                                refine_network = TRUE,
                                                n_neighbors = n_neighbors,
                                                n_components = n_components,
                                                weight_quantile = weight_quantile,
                                                weights_normalization_coef = weights_normalization_coef)

      partial_network <- partial_umap_obj$refined_network

      partial_edges <- partial_umap_obj$umap_out$refined_network$edges

      partial_distances <- partial_umap_obj$umap_out$refined_network$distances
      partial_indexes <- partial_umap_obj$umap_out$refined_network$indexes

      f_ind <- partial_indexes[c,seq(2,min(n_neighbors, dim(f_ind)[2]))]
      f_distances <- partial_distances[c,seq(2,min(n_neighbors, dim(f_ind)[2]))]
      neighbors_c <- rownames(partial_indexes)[f_ind]

    } else {
      f_ind <- f_ind[seq(1,min(n_neighbors, length(f_ind)))]
      f_distances <- f_distances[seq(1,min(n_neighbors, length(f_ind)))]
      neighbors_c <- rownames(indexes)[f_ind]
    }

    xs <- as.numeric(layout[neighbors_c,1])
    ys <- as.numeric(layout[neighbors_c,2])

    if(n_components > 2) {
      zs <- as.numeric(layout[neighbors_c,3])
    }

    weights_to_use <- 1/f_distances
    weights_to_use <- (weights_to_use-min(weights_to_use))/(max(weights_to_use)-min(weights_to_use))

    new_x <- stats::weighted.mean(xs, weights_to_use)
    new_y <- stats::weighted.mean(ys, weights_to_use)

    if(n_components > 2) {
      new_z <- stats::weighted.mean(zs, weights_to_use)
    }

    neighbors_x <- c(neighbors_x, new_x)
    neighbors_y <- c(neighbors_y, new_y)

    if(n_components > 2) {
      neighbors_z <- c(neighbors_z, new_z)
    }

  }

  names(neighbors_x) <- new_clusters
  names(neighbors_y) <- new_clusters
  if(n_components > 2) {
    names(neighbors_z) <- new_clusters
  }

  if(n_components > 2) {
    layout <- layout[,c(1,2,3)]
    new_layout <- rbind(layout, cbind(neighbors_x, neighbors_y, neighbors_z))
  } else {
    layout <- layout[,c(1,2)]
    new_layout <- rbind(layout, cbind(neighbors_x, neighbors_y))
  }

  new_layout <- apply(new_layout, 2, function(i) {
    x <- as.numeric(i)
    names(x) <- rownames(new_layout)
    return(x)
  })

  new_edges_to_use <- edges

  for(c in new_clusters) {

    add_edges <- new_edges[which(new_edges$from == c | new_edges$to == c),]
    add_edges$from.x <- as.numeric(new_layout[add_edges$from, 1])
    add_edges$from.y <- as.numeric(new_layout[add_edges$from, 2])
    add_edges$to.x <- as.numeric(new_layout[add_edges$to, 1])
    add_edges$to.y <- as.numeric(new_layout[add_edges$to, 2])

    new_edges_to_use <- rbind(new_edges_to_use, add_edges)

  }

  new_umap_obj$umap_out$old_layout <- new_umap_obj$umap_out$layout
  new_umap_obj$umap_out$layout <- new_layout

  new_umap_obj$umap_out$refined_network$old_edges <- new_umap_obj$umap_out$refined_network$edges
  new_umap_obj$umap_out$refined_network$edges <- new_edges_to_use

  if(n_components == 2) {
    layout <- cbind(layout, 'Size' = max(2,100/dim(layout)[1]))
  } else {
    layout <- cbind(layout, 'Size' = max(200,10000/dim(layout)[1]))
  }

  layout <- cbind(layout, 'Colors' = col_vector[old_color_attr])

  layout <- cbind(layout, 'Group' = old_color_attr)

  if(!is.null(label_attr)) {
    layout <- cbind(layout, 'SubGroup' = old_label_attr)
  }

  layout <- layout[order(match(layout[,'Group'], names(col_vector)[which(!names(col_vector) %in% new_name)])),]

  if(n_components == 2) {
    if(!is.null(label_attr) && !only_new_points) {
      p <- umap_pointsize(layout = layout,
                          color_attr = layout[,'Group'],
                          label_attr = layout[,'SubGroup'],
                          edges = edges,
                          title = title,
                          new_points_col = new_points_col,
                          legend = legend,
                          legend_name = legend_name,
                          only_new_points = only_new_points,
                          show_edges = show_edges,
                          n_increase = n_increase,
                          alpha = alpha)
    } else {
      p <- umap_pointsize(layout = layout,
                          color_attr = layout[,'Group'],
                          label_attr = NULL,
                          edges = edges,
                          title = title,
                          new_points_col = new_points_col,
                          legend = legend,
                          legend_name = legend_name,
                          only_new_points = only_new_points,
                          show_edges = show_edges,
                          n_increase = n_increase,
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

  if(n_components == 2) {
    new_layout <- cbind(new_layout, 'Size' = max(2,100/dim(new_layout)[1]))
    new_layout[which(color_attr %in% new_name), 'Size'] <- max(2,100/dim(new_layout)[1]) * 2
  } else {
    new_layout <- cbind(new_layout, 'Size' = max(200,10000/dim(new_layout)[1]))
    new_layout[which(color_attr %in% new_name), 'Size'] <- max(200,10000/dim(new_layout)[1]) * 6
  }

  new_layout <- cbind(new_layout, 'Colors' = col_vector[color_attr])

  new_layout <- cbind(new_layout, 'Group' = color_attr)

  if(!is.null(label_attr)) {
    new_layout <- cbind(new_layout, 'SubGroup' = label_attr)
  }

  new_layout <- new_layout[order(match(new_layout[,'Group'], names(col_vector))),]

  if(n_components == 2) {
    if(!is.null(label_attr)) {
      p_new <- umap_pointsize(layout = new_layout,
                              color_attr = new_layout[,'Group'],
                              label_attr = new_layout[,'SubGroup'],
                              edges = new_edges_to_use,
                              title = title,
                              new_points_col = new_points_col,
                              legend = legend,
                              legend_name = legend_name,
                              only_new_points = only_new_points,
                              show_edges = show_edges,
                              n_increase = n_increase,
                              alpha = alpha)
    } else {
      p_new <- umap_pointsize(layout = new_layout,
                              color_attr = new_layout[,'Group'],
                              edges = new_edges_to_use,
                              label_attr = NULL,
                              title = title,
                              new_points_col = new_points_col,
                              legend = legend,
                              legend_name = legend_name,
                              only_new_points = only_new_points,
                              show_edges = show_edges,
                              n_increase = n_increase,
                              alpha = alpha)
    }
  } else {
    p_new <- plotly::plot_ly(x = new_layout[,1],
                             y = new_layout[,2],
                             z = new_layout[,3],
                             type="scatter3d",
                             mode="markers",
                             colors = col_vector,
                             color = new_layout[,'Group'],
                             size = I(as.numeric(new_layout[,'Size'])))

    p_new$x$layoutAttrs <- list('Params' = list('title' = title,
                                                'scene' = list('xaxis' = list('title' = 'UMAP-1'),
                                                               'yaxis' = list('title' = 'UMAP-2'),
                                                               'zaxis' = list('title' = 'UMAP-3'))))

  }

  out <- S4Vectors::List('Original' = list(umap_obj = umap_obj,
                                           umap_plot = p),
                         'New' = list(umap_obj = new_umap_obj,
                                      umap_plot = p_new))

  if(mapping_quality) {
    mapping_quality <- getMappingConfidence(mapped_obj = out,
                                            signatures_cor = signatures_cor)

    out[['MappingQuality']] <- mapping_quality
  }

  if(annotate) {
    new_best_annotation <- nearest_neighbor_annotation(reference_df = reference_df,
                                                       new_clusters = new_clusters,
                                                       umap_obj = out$New,
                                                       color_attr = color_attr[which(!color_attr %in% new_clusters)],
                                                       n_nearest = n_neighbors,
                                                       col_vector = col_vector,
                                                       title = paste0(n_neighbors, ' NearestNeighbors'), ...)

    out[['NearestNeighborsAnnotation']] <- new_best_annotation
  }

  return(out)

}
