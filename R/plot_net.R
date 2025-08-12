#' Create a ggplot of the network
#'
#' @param layout A data frame with the layout informations to plot
#' @param color_attr The color attribute
#' @param edges The edges to visualize
#' @param label_attr A label vector for the clusters
#' @param title A title
#' @param new_points_col The color of newly mapped points
#' @param legend A boolean, if TRUE the legend is shown. Defaults to FALSE
#' @param legend_name The legend name
#' @param only_new_points A boolean, if TRUE only new points are labeled. Defaults
#' to FALSE
#' @param alpha A vector with the alpha information. Defaults to all 1
#' @param show_edges A boolean, if FALSE no edges are shown
#' @param n_increase The increase for newly mapped points. Defaults to 1
#' @param stroke The stroke around the points, defaults to 0.2
#'
#' @return A ggplot
#'
#' @examples
plot_net <- function(layout,
                     color_attr,
                     edges = NULL,
                     label_attr = NULL,
                     title = NULL,
                     new_points_col = "#FF0000",
                     legend = FALSE,
                     legend_name = NULL,
                     only_new_points = FALSE,
                     alpha = NULL,
                     show_edges = TRUE,
                     n_increase = 1,
                     stroke = 0.2) {

  new_edges_col <- new_points_col[1]
  x_shift <- abs(min(as.numeric(layout[,1])))
  y_shift <- abs(min(as.numeric(layout[,2])))
  layout[,1] <- as.numeric(layout[,1]) + x_shift
  layout[,2] <- as.numeric(layout[,2]) + y_shift
  layout <- as.data.frame(layout)
  idxs_new_points <- gtools::mixedorder(unique(layout[,5])[which(unique(layout[,5]) %in% layout[which(layout$Colors %in% new_points_col),5])])
  layout$Colors <- factor(layout[,4],
                          levels = c(unique(layout[,4])[which(!unique(layout[,5]) %in% layout[which(layout$Colors %in% new_points_col), 5])], unique(layout[,4])[which(unique(layout[,4]) %in% layout[which(layout$Colors %in% new_points_col),4])][idxs_new_points]))
  layout$Group <- factor(layout[,5],
                         levels = c(gtools::mixedsort(unique(layout[,5])[which(!unique(layout[,5]) %in% layout[which(layout$Colors %in% new_points_col),5])]), unique(layout[,5])[which(unique(layout[,5]) %in% layout[which(layout$Colors %in% new_points_col),5])][idxs_new_points]))
  if(!is.null(edges)) {
    edges$from.x <- edges$from.x + x_shift
    edges$to.x <- edges$to.x + x_shift
    edges$from.y <- edges$from.y + y_shift
    edges$to.y <- edges$to.y + y_shift
    edges_color <- rep("grey", length(edges$from))
    idxs <- which(edges$from %in% unique(layout[,5])[which(unique(layout[,5]) %in% layout[which(layout$Colors %in% new_points_col),5])])
    idxs <- c(idxs, which(edges$to %in% unique(layout[,5])[which(unique(layout[,5]) %in% layout[which(layout$Colors %in% new_points_col),5])]))
    edges_color[idxs] <- new_edges_col
    edges$weight[idxs] <- edges$weight[idxs] * n_increase
  }
  if(!is.null(label_attr)) {
    layout$SubGroup <- factor(layout[,6],
                              levels = c(gtools::mixedsort(unique(layout[,6])[which(!unique(layout[,6]) %in% layout[which(layout$Colors %in% new_points_col),6])]), gtools::mixedsort(unique(layout[,6])[which(unique(layout[,6]) %in% layout[which(layout$Colors %in% new_points_col),6])])))
  }
  if(is.null(alpha)) {
    alpha <- rep(1, dim(layout)[1])
  }
  layout$Alpha <- alpha
  colors_to_use <- unique(as.vector(layout[,4]))
  if(any(colors_to_use %in% new_points_col)) {
    colors_to_use <- colors_to_use[-which(colors_to_use %in% new_points_col)]
    if(length(new_points_col) == 1) {
      colors_to_use <- c(colors_to_use, rep(new_points_col,
                                            dim(layout[which(layout$Colors == new_points_col),])[1]))
    }
    else {
      colors_to_use <- c(colors_to_use, new_points_col[as.vector(unique(layout[,5])[which(unique(layout[,5]) %in% names(new_points_col))])])
    }
  }
  names(colors_to_use) <- unique(layout[,5])
  if(is.null(names(new_points_col)) && !is.null(new_points_col)) {
    if(length(new_points_col) == 1) {
      new_points_col <- rep(new_points_col, length(which(layout[,4] %in% new_points_col)))
      names(new_points_col) <- layout[which(layout[,4] %in%
                                              new_points_col), 5]
    } else {
      names(new_points_col) <- layout[match(new_points_col, layout[,4]),5]
    }
  }
  if(any(!is_number(as.vector(layout[,5])))) {
    continuous <- FALSE
  } else {
    continuous <- TRUE
  }
  if(continuous) {
    to_use <- rep(NA, length(as.vector(layout[,5])))
    names(to_use) <- as.vector(layout[,5])
    to_use[which(is_number(as.vector(layout[,5])))] <- as.numeric(as.vector(layout[,5]))
    layout[which(layout[,4] %in% new_points_col),4] <- sort(as.vector(layout[which(layout[,4] %in% new_points_col),4]))
    if(length(unique(new_points_col)) > 1) {
      names(new_points_col) <- layout[match(new_points_col,layout[,4]), ]
    }
    p <- ggplot2::ggplot()
    if(!is.null(edges) & show_edges) {
      p <- p + ggplot2::geom_segment(data = edges, ggplot2::aes(x = from.x,
                                                                xend = to.x,
                                                                y = from.y,
                                                                yend = to.y),
                                     colour = edges_color,
                                     na.rm = TRUE,
                                     linewidth = edges$weight*(min(1, stroke+0.5)),
                                     alpha = edges$weight/max(edges$weight))
    }
    p <- p + ggplot2::geom_point(data = (layout),
                                 mapping = ggplot2::aes(x = as.numeric((layout)[,1]),
                                                        y = as.numeric((layout)[,2]),
                                                        fill = to_use,
                                                        stroke = stroke*(alpha**2)),
                                 alpha = alpha,
                                 shape = 21,
                                 size = as.numeric((layout)[,3]),
                                 show.legend = legend) +
      ggplot2::scale_colour_gradientn(colours = colors_to_use[which(!colors_to_use %in% new_points_col)],
                                      name = legend_name,
                                      na.value = new_points_col) +
      ggplot2::scale_fill_gradientn(colours = colors_to_use[which(!colors_to_use %in% new_points_col)],
                                    name = legend_name,
                                    na.value = new_points_col) +
      ggplot2::theme(axis.text = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     legend.key = ggplot2::element_blank(),
                     legend.position = "right") +
      ggplot2::scale_alpha(guide = "none")
  }
  else {
    p <- ggplot2::ggplot()
    if(!is.null(edges) & show_edges) {
      p <- p + ggplot2::geom_segment(data = edges,
                                     ggplot2::aes(x = from.x,
                                                  xend = to.x,
                                                  y = from.y,
                                                  yend = to.y),
                                     colour = edges_color,
                                     na.rm = TRUE,
                                     linewidth = edges$weight*(min(1, stroke+0.5)),
                                     alpha = edges$weight/max(edges$weight))
    }
    p <- p + ggplot2::geom_point(data = (layout),
                                 mapping = ggplot2::aes(x = as.numeric((layout)[,1]),
                                                        y = as.numeric((layout)[,2]),
                                                        color = layout[,5],
                                                        stroke = stroke*(layout[,dim(layout)[2]])**2,
                                                        fill = layout[,5]), alpha = layout[,dim(layout)[2]],
                                 shape = 21,
                                 size = as.numeric((layout)[,3]),
                                 show.legend = legend) +
      ggplot2::scale_color_manual(values = rep("black", length(colors_to_use)),
                                  guide = ggplot2::guide_legend(title.position = "top",
                                                                title.hjust = 0.5,
                                                                title = legend_name)) +
      ggplot2::scale_fill_manual(values = colors_to_use,
                                 guide = ggplot2::guide_legend(title.position = "top",
                                                               title.hjust = 0.5,
                                                               title = legend_name,
                                                               override.aes = list(size = 5))) +
      ggplot2::theme(axis.text = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     legend.key = ggplot2::element_blank(),
                     legend.position = "right") +
      ggplot2::scale_alpha(guide = "none")
  }
  if(!is.null(title)) {
    p <- p + ggplot2::ggtitle(title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
                                                        size = max(15, min(30, dim(layout)[1]/20)),
                                                        face = "bold"))
  }
  if(is.null(label_attr)) {
    return(p)
  }
  if(all(color_attr == label_attr)) {
    all_cx <- c()
    all_cy <- c()
    all_colors <- c()
    for (g in unique(color_attr)) {
      c <- layout[which(color_attr == g), ]
      if(continuous && !any((as.vector(c[,4])) %in%
                            new_points_col)) {
        mean_all <- mean(as.numeric(as.vector(c[,5])))
        all_days <- as.numeric(as.vector(c[,5]))
        distances <- abs(all_days - mean_all)
        c <- c[which(as.numeric(as.vector(c[,5])) == all_days[which.min(distances)]),]
      }
      if(dim(c)[1] == 1) {
        all_cx <- c(all_cx, as.numeric(c[1]))
        all_cy <- c(all_cy, as.numeric(c[2]))
        all_colors <- c(all_colors, as.vector(c[[4]]))
      }
      else {
        v_x <- as.numeric(c[,1])
        v_y <- as.numeric(c[,2])
        mx <- stats::median(v_x)
        my <- stats::median(as.numeric(c[,2]))
        external_point <- c(mx, my)
        points <- apply(c[,c(1, 2)], 2, as.numeric)
        dif_x <- (points[,1] - external_point[1])^2
        dif_y <- (points[,2] - external_point[2])^2
        distances <- sqrt(rowSums(cbind(dif_x, dif_y)))
        closest_point_index <- which.min(distances)
        closest_point <- c[closest_point_index, ]
        all_cx <- c(all_cx, as.numeric(closest_point[1]))
        all_cy <- c(all_cy, as.numeric(closest_point[2]))
        all_colors <- c(all_colors, as.vector(closest_point[[4]]))
      }
    }
    if(any(new_points_col %in% all_colors)) {
      all_colors <- all_colors[which(!all_colors %in%
                                       new_points_col)]
      l_f <- layout[which(layout[,4] %in% new_points_col),]
      if(is.vector(l_f)) {
        l_f_groups <- l_f[5]
      }
      else {
        l_f_groups <- l_f[,5]
      }
      if(length(new_points_col) == 1) {
        all_colors <- c(all_colors, rep(new_points_col, length(unique(l_f_groups))))
      }
      else {
        all_colors <- c(all_colors, new_points_col[unique(label_attr)[which(unique(label_attr) %in% names(new_points_col))]])
      }
    }
    idx <- which(all_colors %in% new_points_col)
    label <- unique(color_attr)
    if(length(idx) > 0) {
      if(continuous && length(new_points_col) > 1) {
        all_colors[idx] <- "black"
      }
      if(only_new_points) {
        all_cx <- all_cx[idx]
        all_cy <- all_cy[idx]
        all_colors <- all_colors[idx]
        label <- label[idx]
      }
    }
    all_colors_labels <- all_colors
    if(any(is.na(all_colors_labels))) {
      all_colors_labels[which(is.na(all_colors_labels))] <- "grey90"
    }
    p <- p + ggrepel::geom_label_repel(ggplot2::aes(x = all_cx,
                                                    y = all_cy,
                                                    label = label),
                                       color = all_colors_labels,
                                       box.padding = 0.5,
                                       max.overlaps = Inf,
                                       show.legend = FALSE)
  }
  else {
    all_cx <- c()
    all_cy <- c()
    all_colors <- c()
    for (l in unique(label_attr)) {
      c <- layout[which(label_attr == l), ]
      if(continuous && !any((as.vector(c[,4])) %in% new_points_col)) {
        mean_all <- mean(as.numeric(as.vector(c[,5])))
        all_days <- as.numeric(as.vector(c[,5]))
        distances <- abs(all_days - mean_all)
        c <- c[which(as.numeric(as.vector(c[,5])) == all_days[which.min(distances)]),]
      }
      if(dim(c)[1] == 1) {
        all_cx <- c(all_cx, as.numeric(c[1]))
        all_cy <- c(all_cy, as.numeric(c[2]))
        all_colors <- c(all_colors, as.vector(c[[4]]))
      }
      else {
        v_x <- as.numeric(c[,1])
        v_y <- as.numeric(c[,2])
        mx <- stats::median(v_x)
        my <- stats::median(as.numeric(c[,2]))
        external_point <- c(mx, my)
        points <- apply(c[,c(1, 2)], 2, as.numeric)
        dif_x <- (points[,1] - external_point[1])^2
        dif_y <- (points[,2] - external_point[2])^2
        distances <- sqrt(rowSums(cbind(dif_x, dif_y)))
        closest_point_index <- which.min(distances)
        closest_point <- c[closest_point_index, ]
        all_cx <- c(all_cx, as.numeric(closest_point[1]))
        all_cy <- c(all_cy, as.numeric(closest_point[2]))
        all_colors <- c(all_colors, as.vector(closest_point[[4]]))
      }
    }
    if(any(new_points_col %in% all_colors)) {
      all_colors <- all_colors[which(!all_colors %in%
                                       new_points_col)]
      l_f <- layout[which(layout[,4] %in% new_points_col),]
      l_f_groups <- ifelse(is.vector(lf), l_f[5], l_f[,5])
      if(length(new_points_col) == 1) {
        all_colors <- c(all_colors, rep(new_points_col, length(unique(l_f_groups))))
      }
      else {
        all_colors <- c(all_colors, new_points_col[unique(label_attr)[which(unique(label_attr) %in% names(new_points_col))]])
      }
    }
    idx <- which(all_colors %in% new_points_col)
    label <- unique(label_attr)
    if(length(idx) > 0) {
      if(continuous && length(new_points_col) > 1) {
        all_colors[idx] <- "black"
      }
      if(only_new_points) {
        all_cx <- all_cx[idx]
        all_cy <- all_cy[idx]
        all_colors <- all_colors[idx]
        label <- label[idx]
      }
    }
    all_colors_labels <- all_colors
    if(any(is.na(all_colors_labels))) {
      all_colors_labels[which(is.na(all_colors_labels))] <- "grey90"
    }
    p <- p + ggrepel::geom_label_repel(ggplot2::aes(x = all_cx,
                                                    y = all_cy,
                                                    label = label),
                                       color = all_colors_labels,
                                       box.padding = 0.5,
                                       max.overlaps = Inf,
                                       show.legend = FALSE)
  }
  return(p)
}
