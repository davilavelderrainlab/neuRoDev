#' Plot a score onto the network
#'
#' @inheritParams plotNetwork
#' @param genes The genes to consider for the score. If multiple genes are given,
#' the score will be their average.
#' @param score The score to plot. If NULL (default), the logcounts of the
#' pseudobulks are used
#' @param smooth A boolean, if TRUE (default), the scores are smoothed to consider
#' the network structure
#' @param n_nearest The number of nearest neighbors to consider for smoothing.
#' Defaults to 15
#' @param normalize A boolean, if TRUE (default) and score is NULL, it normalizes
#' the pseudobulk expression computing a row-wise z-score
#' @param na.vec The vector to define if certain clusters are not to be shown
#' @param palette The color palette to use. Defaults to `blues9`
#' @param title The title of the plot
#' @param show_edges A boolean, if FALSE no edges are shown. Defaults to TRUE
#' @param stroke The stroke of the points, defaults to 0.5
#' @param fix_alpha A boolean, if TRUE all points will have alpha 1, otherwise
#' it will be proportional to the score. Defaults to FALSE
#'
#' @return A ggplot of the network colored based on the score
#' @export
#'
#' @examples
#' m <- matrix(sample(seq(1,10, length.out=10000), 15000*100, replace = TRUE), ncol = 100)
#' rownames(m) <- paste0('Gene-', seq(1,15000))
#' colnames(m) <- paste0('Col-', seq(1,100))
#' net <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = m))
#' net$SubClass <- rep(c('A', 'B', 'C', 'D'), each = 25)
#' subclass_palette <- c('A' = 'red', 'B' = 'blue', 'C' = 'green', 'D' = 'yellow')
#' net$SubClass_colors <- subclass_palette[net$SubClass]
#' net$X_coord <- sample(seq(1,2, length.out = 1000), size = ncol(net), replace = TRUE)
#' net$Y_coord <- sample(seq(1,2, length.out = 1000), size = ncol(net), replace = TRUE)
#' edges_from <- sample(colnames(net), size = 200, replace = TRUE)
#' edges_to <- sample(colnames(net), size = 200, replace = TRUE)
#' edges_from_x <- net$X_coord[match(edges_from, colnames(net))]
#' edges_from_y <- net$Y_coord[match(edges_from, colnames(net))]
#' edges_to_x <- net$X_coord[match(edges_to, colnames(net))]
#' edges_to_y <- net$Y_coord[match(edges_to, colnames(net))]
#' edges_weight <- sample(seq(0,1, length.out=1000), length(edges_from), replace = TRUE)
#' edges_df <- data.frame('from' = edges_from, 'to' = edges_to, 'weight' = edges_weight,
#' 'from.x' = edges_from_x,
#' 'from.y' = edges_from_y,
#' 'to.x' = edges_to_x,
#' 'to.y' = edges_to_y)
#' net@metadata$network$edges <- edges_df
#' plotNetworkScore(net, 'Gene-1')
plotNetworkScore <- function(net,
                             genes = NULL,
                             score = NULL,
                             label_attr = NULL,
                             smooth = FALSE,
                             n_nearest = 15,
                             normalize = TRUE,
                             na.vec = NULL,
                             palette = grDevices::blues9,
                             title=NULL,
                             show_edges = TRUE,
                             stroke = 0.5,
                             fix_alpha = FALSE,
                             max_size = 3) {

  if(is.character(label_attr) && length(label_attr) == 1 && label_attr %in% colnames(SingleCellExperiment::colData(net))) {
    label_attr <- SingleCellExperiment::colData(net)[,label_attr]
  }

  label_attr <- as.vector(label_attr)

  if(is.null(genes)) {
    genes <- rownames(net)
  } else {
    if(!all(genes %in% rownames(net))) {
      message(paste0(genes[which(!genes %in% rownames(net))], collapse = '-'), ' not in `net`')
      genes <- genes[which(genes %in% rownames(net))]
      if(length(genes) == 0) {
        return('No gene in `net`')
      }
    }
    genes <- genes[which(genes %in% rownames(net))]
  }

  if(is.null(score)) {
    exp_mat <- SingleCellExperiment::logcounts(net)
    if(normalize) {
      exp_mat <- t(apply(exp_mat, 1, function(v) {(v-mean(v))/stats::sd(v)}))
    }
    score <- exp_mat[genes,]
  } else {
    if(is.matrix(score)) {
      if(!all(genes %in% rownames(score))) {
        message(paste0(genes[which(!genes %in% rownames(score))], collapse = '-'), ' not in `score`')
        genes <- genes[which(genes %in% rownames(score))]
        if(length(genes) == 0) {
          return('No gene in `score`')
        }
      }
      if(normalize) {
        score <- t(apply(score, 1, function(v) {(v-mean(v))/stats::sd(v)}))
      }
      score <- score[genes,]
    }
  }

  if(!is.vector(score)) {
    score <- Matrix::colMeans(score)
  }

  if(is.null(names(score))) {
    names(score) <- colnames(net)
  } else {
    score <- score[colnames(net)]
  }

  edges <- net@metadata$network$edges

  if(smooth) {
    network <- igraph::graph_from_adjacency_matrix(as.matrix(net@metadata$network$adj_matrix),
                                                   weighted = TRUE)
    new_score <- network_smoothing(network = network,
                                   scores = score[igraph::V(network)$name],
                                   n_nearest_smooth = n_nearest)

    new_score <- new_score[names(score)]
    score <- new_score

  }

  if(is.null(na.vec)) {
    na.vec <- rep(1, length(score))
  }

  breakpoints <- stats::quantile(seq((min(score[which(!is.na(na.vec))]) * rounder(1/min(score[which(!is.na(na.vec))]))),
                                     (max(score[which(!is.na(na.vec))]) * rounder(1/min(score[which(!is.na(na.vec))])))),
                                 probs = seq(0, 1, length.out = length(palette)+1))
  breakpoints[which.min(breakpoints)] <- breakpoints[which.min(breakpoints)] - 1
  breakpoints[which.max(breakpoints)] <- breakpoints[which.max(breakpoints)] + 1
  col_vector <- as.vector(cut(score * rounder(1/min(score)), breaks = breakpoints, labels = palette))

  if(any(is.na(na.vec))) {
    col_vector[which(is.na(na.vec))] <- NA
  }
  names(col_vector) <- score
  col_vector <- col_vector[order(score)]

  color_attr <- label_attr

  layout <- cbind(net$X_coord, net$Y_coord)
  rownames(layout) <- colnames(net)
  size <- rep(max(max_size, 100/dim(layout)[1]), dim(layout)[1])
  layout <- cbind(layout, Size = size)
  layout <- cbind(layout, Colors = col_vector[as.character(score)])
  layout <- cbind(layout, Group = col_vector[as.character(score)])
  if(!is.null(label_attr)) {
    layout <- cbind(layout, SubGroup = color_attr)
  }

  alpha <- min_max_normalize(score)

  if(any(is.na(col_vector))) {
    alpha[order(score)][which(is.na(col_vector))] <- 0
  }
  min_maxed_score <- alpha
  if(fix_alpha) {
    alpha <- NULL
  }

  if(is.null(title)) {
    if(length(genes) != nrow(net)) {
      title <- paste(genes[seq(1,min(length(genes), 5))], collapse = "-")
    }
  }

  weights <- apply(edges, 1, function(i) {
    i[["weight"]] <- as.numeric(i[["weight"]]) * min(c(as.numeric(min_maxed_score[i[["from"]]]),
                                                       as.numeric(min_maxed_score[i[["to"]]])))
    return(i[["weight"]])
  })
  edges$weight <- as.numeric(weights)

  if(!is.null(label_attr)) {
    p <- plot_net(layout = layout,
                  color_attr = layout[,"Colors"],
                  label_attr = layout[,"SubGroup"],
                  title = title,
                  legend = FALSE,
                  only_new_points = FALSE,
                  alpha = alpha,
                  edges = edges,
                  show_edges = show_edges,
                  stroke = stroke)
  } else {
    p <- plot_net(layout = layout,
                  color_attr = layout[,"Colors"],
                  label_attr = label_attr,
                  title = title,
                  legend = FALSE,
                  only_new_points = FALSE,
                  alpha = alpha,
                  edges = edges,
                  show_edges = show_edges,
                  stroke = stroke)
  }

  return(p)

}
