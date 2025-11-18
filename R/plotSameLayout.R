#' Plot a mapped dataset conserving the same layout as the original reference
#' network
#'
#' @inheritParams mapNetwork
#' @param new_cor The correlation matrix between the profiles and the reference
#' pseudobulks
#' @param no_label A boolean, if TRUE the clusters' labels will not be shown.
#' Defaults to FALSE
#' @param new_points_col The color of the mapped points
#' @param max_size The maximum size of the dots in the network visualization
#' @param only_new_points A boolean, if TRUE only the labels of the new points
#' are shown. Defaults to FALSE
#' @param annotate A boolean, if TRUE (default), it returns also the
#' annotation scores for the mapped points
#' @param order_names A boolean, if TRUE the names in annotateMapping
#' barplot are ordered. Defaults to FALSE
#' @param plot A boolean, if TRUE the mapping plot is produce, otherwise only the
#' mapping annotation is given. Defaults to TRUE
#' @param ... Additional parameters for `plot_net`
#'
#' @return An `S4Vectors` list containing the correlation matrix between the
#' profiles and the reference pseudobulks (`new_cor`), a data.frame containing
#' the edges information to connect the new profiles with the pseudobulks in the
#' network (`new_edges`), a matrix containing the layout information, including
#' the x and y coordinates (`new_layout`), the new plot (`new_plot`), a mapping
#' quality assessment (if mapping_quality = TRUE) (`mapping_quality`), and
#' the annotation scores (if annotate = TRUE) (`annotation`)
#'
#' @examples
#' m <- matrix(sample(seq(1,10, length.out=10000), 15000*100, replace = TRUE), ncol = 100)
#' rownames(m) <- paste0('Gene-', seq(1,15000))
#' colnames(m) <- paste0('Col-', seq(1,100))
#' net <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = m))
#' net$SubClass <- rep(c('A', 'B', 'C', 'D'), each = 25)
#' subclass_palette <- c('A' = 'red', 'B' = 'blue', 'C' = 'green', 'D' = 'yellow')
#' net$SubClass_color <- subclass_palette[net$SubClass]
#' net$Stages <- rep(c('S1', 'S2', 'S3', 'S4'), each = 25)
#' stages_palette <- c('S1' = 'pink', 'S2' = 'orange', 'S3' = 'violet', 'S4' = 'black')
#' net$Stages_color <- stages_palette[net$Stages]
#' net$X_coord <- sample(seq(1,2, length.out = 1000), size = ncol(net), replace = TRUE)
#' net$Y_coord <- sample(seq(1,2, length.out = 1000), size = ncol(net), replace = TRUE)
#' edges_from <- sample(colnames(net), size = 200, replace = TRUE)
#' edges_to <- sample(colnames(net), size = 200, replace = TRUE)
#' edges_from_x <- net$X_coord[match(edges_from, colnames(net))]
#' edges_from_y <- net$Y_coord[match(edges_from, colnames(net))]
#' edges_to_x <- net$X_coord[match(edges_to, colnames(net))]
#' edges_to_y <- net$Y_coord[match(edges_to, colnames(net))]
#' edges_weight <- sample(seq(0,1, length.out=1000), length(edges_from), replace = TRUE)
#' edges_df <- data.frame('from' = edges_from,
#' 'to' = edges_to,
#' 'weight' = edges_weight,
#' 'from.x' = edges_from_x,
#' 'from.y' = edges_from_y,
#' 'to.x' = edges_to_x,
#' 'to.y' = edges_to_y)
#' net@metadata$network$edges <- edges_df
#' SummarizedExperiment::rowData(net)$informative <- sample(c(TRUE, FALSE), size = nrow(net),
#' replace = TRUE)
#' new_profiles <- matrix(sample(seq(1,10, length.out=10000), nrow(net)*10, replace = TRUE), ncol = 10)
#' rownames(new_profiles) <- rownames(net)
#' colnames(new_profiles) <- paste0('NewCol-', seq(1,10))
#' common_genes <- intersect(rownames(net)[SingleCellExperiment::rowData(net)$informative],
#' rownames(new_profiles))
#' new_cor <- stats::cor(t(apply(as.matrix(SingleCellExperiment::logcounts(net)[common_genes,]),
#' 1, function(v) {(v-mean(v))/stats::sd(v)})),new_profiles[common_genes,])
#' neuRoDev:::plotSameLayout(net, new_cor)
plotSameLayout <- function(net,
                           new_cor,
                           color_attr = 'SubClass',
                           label_attr = 'SubClass',
                           col_vector = NULL,
                           no_label = FALSE,
                           new_points_col = "#FF0000",
                           max_size = 3,
                           only_new_points = FALSE,
                           annotate = TRUE,
                           n_nearest = 15,
                           new_name = NULL,
                           order_names = FALSE,
                           plot = TRUE,
                           ...) {

  if(is.null(col_vector)) {
    col_vector <- paste0(color_attr, '_color')
  }

  orig_col_attr <- color_attr
  orig_label_attr <- label_attr
  orig_col_vec <- col_vector

  if(!plot | ncol(new_cor) > 1000) {
    new_best_annotation <- annotateMapping(net = net,
                                           new_cor = new_cor,
                                           color_attr = orig_col_attr,
                                           col_vector = orig_col_vec,
                                           n_nearest = n_nearest,
                                           order_names = order_names)

    out <- S4Vectors::List(new_cor = new_cor,
                                annotation = new_best_annotation)

    return(out)
  }

  if(is.character(label_attr) && length(label_attr) == 1 && label_attr %in% colnames(SingleCellExperiment::colData(net))) {
    label_attr <- SingleCellExperiment::colData(net)[,label_attr]
  }
  if(is.character(color_attr) && length(color_attr) == 1 && color_attr %in% colnames(SingleCellExperiment::colData(net))) {
    color_attr <- SingleCellExperiment::colData(net)[,color_attr]
  }
  if(is.character(col_vector) && length(col_vector) == 1 && col_vector %in% colnames(SingleCellExperiment::colData(net))) {
    col_vector <- SingleCellExperiment::colData(net)[,col_vector]
  }
  label_attr <- as.vector(label_attr)
  color_attr <- as.vector(color_attr)
  names(col_vector) <- color_attr

  new_clusters <- colnames(new_cor)
  if(is.null(new_name)) {
    new_name <- new_clusters
  }

  if(is.null(label_attr) & !no_label) {
    label_attr <- color_attr
  }
  if(no_label) {
    label_attr <- NULL
  }

  if(is.null(names(new_points_col))) {
    if(length(new_points_col) == 1) {
      new_points_col <- rep(new_points_col,
                            length(new_clusters)+length(new_name))
    }
    names(new_points_col) <- c(new_clusters, new_name)
  }

  if(length(color_attr) == ncol(net)) {
    old_color_attr <- color_attr
    if(length(new_name) == 1) {
      color_attr <- c(color_attr, rep(new_name, nrow(new_cor)))
    } else {
      color_attr <- c(color_attr, new_name)
    }
  } else {
    old_color_attr <- color_attr[which(!color_attr %in% new_name)]
  }

  if(!is.null(label_attr)) {
    if(length(label_attr) == ncol(net)) {
      old_label_attr <- label_attr
      if(length(new_name) == 1) {
        label_attr <- c(label_attr, rep(new_name, dim(new_cor)[1]))
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
    new_clusters[which(new_clusters %in% c(old_label_attr,old_color_attr))] <- paste0("New-", new_clusters[which(new_clusters %in% c(old_label_attr, old_color_attr))])
  }
  colnames(new_cor) <- new_clusters

  layout <- cbind(net$X_coord, net$Y_coord)
  rownames(layout) <- colnames(net)
  edges <- net@metadata$network$edges
  network <- net@metadata$network$adj_matrix

  if(length(unique(col_vector)) == length(unique(old_color_attr))) {
    if(!is.null(names(col_vector))) {
      old_names <- names(col_vector)
      col_vector <- c(col_vector, new_points_col[c(new_clusters, new_name)])
      names(col_vector) <- c(old_names, new_clusters, new_name)
    } else {
      col_vector <- c(col_vector, new_points_col[c(new_clusters, new_name)])
      names(col_vector) <- c(unique(old_color_attr), new_clusters, new_name)
    }
  }

  center <- c(mean(as.numeric(layout[,1])), mean(as.numeric(layout[,2])))

  annot <- rownames(new_cor)
  layout_x <- as.numeric(layout[,1])
  layout_y <- as.numeric(layout[,2])

  neighbors_coords <- do.call(rbind,lapply(seq_len(ncol(new_cor)), function(j) {
      col <- new_cor[, j]

      ord <- order(col, decreasing = TRUE)
      idx <- ord[seq_len(n_nearest)]

      xs <- layout_x[idx]
      ys <- layout_y[idx]

      sim <- col[idx]

      w <- sim / sum(sim)
      weighted_x <- sum(xs * w)
      weighted_y <- sum(ys * w)

      avg_sim <- mean(sim)
      factor <- min(1, 2 * avg_sim)

      return((1 - factor) * center + factor * c(weighted_x, weighted_y))
  }))

  rownames(neighbors_coords) <- colnames(new_cor)

  new_layout <- rbind(layout, neighbors_coords)

  new_edges <- do.call(rbind, lapply(colnames(new_cor), function(coln) {
    i <- new_cor[,coln]
    n <- rownames(new_cor)[order(i, decreasing = TRUE)[seq(1,n_nearest)]]
    cors <- i[n]
    edge <- data.frame('from' = rep(coln, n_nearest),
                       'to' = n,
                       'weight' = cors/min(2,(max(1.25, n_nearest/20))),
                       'from.x' = new_layout[coln, 1],
                       'from.y' = new_layout[coln, 2],
                       'to.x' = new_layout[n, 1],
                       'to.y' = new_layout[n, 2])
    return(edge)
  }))

  rownames(new_edges) <- seq(1, nrow(new_edges))

  new_edges_to_use <- rbind(edges, new_edges)

  new_layout <- cbind(new_layout, Size = max(max_size, 100/dim(new_layout)[1]))
  new_layout[which(color_attr %in% new_name), "Size"] <- max(max_size, 100/dim(new_layout)[1])*1.5
  new_layout <- cbind(new_layout, Colors = col_vector[color_attr])
  new_layout <- cbind(new_layout, Group = color_attr)
  if(!is.null(label_attr)) {
    new_layout <- cbind(new_layout, SubGroup = label_attr)
  }
  new_layout <- new_layout[order(match(new_layout[, "Group"],names(col_vector))),]

  if(!is.null(label_attr)) {
    p_new <- plot_net(layout = new_layout,
                      color_attr = new_layout[, "Group"],
                      label_attr = new_layout[, "SubGroup"],
                      edges = new_edges_to_use,
                      new_points_col = new_points_col,
                      only_new_points = only_new_points, ...)
  }
  else {
    p_new <- plot_net(layout = new_layout,
                      color_attr = new_layout[, "Group"],
                      label_attr = NULL,
                      edges = new_edges_to_use,
                      new_points_col = new_points_col,
                      only_new_points = only_new_points, ...)
  }

  out <- S4Vectors::List(new_cor = new_cor,
                              new_edges = new_edges_to_use,
                              new_layout = new_layout,
                              new_plot = p_new)

  if(annotate) {
    new_best_annotation <- annotateMapping(net = net,
                                           new_cor = new_cor,
                                           color_attr = orig_col_attr,
                                           col_vector = orig_col_vec,
                                           n_nearest = n_nearest,
                                           order_names = order_names)
    out[["annotation"]] <- new_best_annotation
  }

  return(out)

}
