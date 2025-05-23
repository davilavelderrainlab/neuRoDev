#' Nearest neighbor annotation
#'
#'  Given an annotated reference dataframe reference_df,
#'  a set of new clusters to analyze in new_clusters, a umap object as run
#'  with add_to_reference, which contains the new clusters to
#'  analyze, and a given annotation label (and possibly a
#'  sub_color_attr), it returns the annotation (and sub_annotation) of
#'  each cluster based on its nearest neighbors. With n_nearest you can set
#'  how many nearest neighbours to consider (it will consider first the
#'  closest ones). You can provide a col_vector and sub_col_vector palettes
#'  for the final barplots (otherwise it will create a palette).
#' @inheritParams umap_plot_same_layout
#' @inheritParams add_to_reference
#' @param col_vector A named color palette vector for the annotations
#' @param n_nearest The number of nearest neighbors to consider. It defaults
#' to `((dim(distance_matrix)[2])-1)`, where distance_matrix is given by:
#' `umap_obj$umap_out$knn$distances` or
#' `umap_obj$umap_out$refined_network$distances`
#' @param title The title to add to the plots
#' @param to_exclude The labels inside `reference_df[,color_attr]`
#' that the user doesn't want to consider.
#' @param compute_means A boolean variable. If TRUE, it computes the mean instead
#' of the sum of values belonging to a given group.
#' @param color_attr The labels to use for the annotation
#'
#' @return A list that contains, for each new cluster, the annotations of the
#' neighbors, the most frequent annotation, the highest value in terms of
#' weighted frequency (based on the distance of the neighbor), a Heatmap and a
#' barplot and the final number of nearest neighbors considered. If a sub
#' annotation is provided, the same will be given for the sub annotation as well
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
#' res <- add_to_reference(annotated_M,
#' new_M,
#' annotated_M$`Best.Assignment`)
#' umap_obj <- res$New
#' new_clusters <- rownames(new_M)
#' nearest_neighbor_annotation(annotated_M,
#' new_clusters,
#' umap_obj,
#' color_attr = 'Best.Assignment')
nearest_neighbor_annotation <- function(reference_df,
                                        new_clusters,
                                        umap_obj,
                                        color_attr,
                                        col_vector=NULL,
                                        n_nearest=NULL,
                                        title=NULL,
                                        to_exclude=NULL,
                                        compute_means=FALSE) {

  if(is.character(color_attr) & length(color_attr) == 1) {
    color_attr <- reference_df[,color_attr]
  }

  idx_annotation <- which(apply(reference_df, 2, function(i) {
    all(i == color_attr)
  }))

  annotation_legend <- colnames(reference_df)[idx_annotation]

  if('umap_obj' %in% names(umap_obj)) {
    umap_obj <- umap_obj$umap_obj
  }

  if(!is.null(umap_obj$refined_network)) {
    index_matrix <- umap_obj$umap_out$refined_network$indexes
    distance_matrix <- umap_obj$umap_out$refined_network$distances
  } else {
    index_matrix <- umap_obj$umap_out$knn$indexes
    distance_matrix <- umap_obj$umap_out$knn$distances
  }

  n_nearest <- min(n_nearest, ((dim(distance_matrix)[2])-1))

  min_n_nearest <- c()
  for(c in new_clusters) {
    if(!c %in% rownames(index_matrix)) {
      c <- paste0('New-', c)
      new_clusters[which(new_clusters == c)] <- paste0('New-', c)
    }
    all_n <- rownames(index_matrix)[index_matrix[c,seq(2,dim(index_matrix)[2])]]
    idxs <- which(!(all_n %in% new_clusters))
    min_n_nearest <- c(min_n_nearest, length(idxs))
  }

  n_nearest <- min(min(min_n_nearest), n_nearest)

  if(n_nearest > 0.9*dim(reference_df)[1] & compute_means == FALSE) {
    message('The number of nearest neighbours
            (computed as min(n_nearest, ((dim(distance_matrix)[2])-1)))
            is higher than 90% of the total clusters. It should be better to
            compute the means instead of the sum (compute_means = TRUE), to
            avoid biases of the base structure of the clusters')
  }

  if(!is.null(to_exclude)) {
    idxs_to_remove <- which(color_attr %in% to_exclude)
  } else {
    idxs_to_remove <- NULL
  }

  annotation_tables <- list()

  clusters_to_exclude <- c(new_clusters, rownames(reference_df)[idxs_to_remove])

  for(c in new_clusters) {
    dist <- distance_matrix[c,seq(2,dim(distance_matrix)[2])]
    all_n <- rownames(index_matrix)[index_matrix[c,seq(2,dim(index_matrix)[2])]]
    idxs <- which(!(all_n %in% clusters_to_exclude))
    all_n <- all_n[idxs]
    dist <- dist[idxs]
    all_n <- all_n[seq(1,n_nearest)]
    dist <- dist[seq(1,n_nearest)]
    reference_df_f <- reference_df[all_n[which(all_n %in% rownames(reference_df))],]

    names(dist) <- reference_df_f[,idx_annotation]

    annotations <- value_table(dist,
                               perc = TRUE,
                               reciprocal = TRUE,
                               compute_means = compute_means)

    annotation_tables[[c]] <- annotations
  }

  max_annotation <- unlist(lapply(annotation_tables, function(i) {
    names(i)[which.max(i)]
  }))
  max_value <- unlist(lapply(annotation_tables, function(i) {
    i[which.max(i)]
  }))
  names(max_value) <- max_annotation

  m <- get_matrix_from_list(annotation_tables)

  if(length(unique(as.vector(m)))==1) {
    return('All clusters have the same scores, you should probably change what
           you map or the label attribute used.')
  }

  h1 <- ComplexHeatmap::Heatmap(m,
                col = grDevices::blues9,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                name = 'Neighborhood score\nAnnotation')

  if(is.null(col_vector)) {
    col_vector <- Polychrome::glasbey.colors(length(colnames(m))+1)
    col_vector <- col_vector[seq(2,length(col_vector))]
    names(col_vector) <- colnames(m)
  }

  df <- reshape2::melt(m)

  df$Var2 <- factor(df$Var2, levels = gtools::mixedsort(unique(as.vector(df$Var2)), decreasing = TRUE))

  plot <- ggplot2::ggplot(df, ggplot2::aes(x = Var1,
                         y = value,
                         fill = Var2)) +
    ggplot2::geom_bar(stat = "identity")    +
    ggplot2::labs(x = "",
         y = "Values",
         fill = annotation_legend,
         title = title) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_manual(values = col_vector[match(df$Var2, names(col_vector))]) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(hjust = 0.5,
                                    size = 20),
          axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))

  return(S4Vectors::SimpleList('Annotations' = annotation_tables,
                               'AnnotationsMatrix' = m,
                               'Best.Annotation' = max_annotation,
                               'Best.Value' = max_value,
                               'Heatmap-Annotation' = h1,
                               'Barplot-Annotation' = plot,
                               'Number_of_NN' = n_nearest))
}
