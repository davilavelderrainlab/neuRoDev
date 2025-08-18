#' Map new profiles to the eTraces
#'
#' @inheritParams plot_eTrace
#' @param mapped_obj The result of `mapNetwork`
#' @param n_nearest The number of nearest neighbors to consider to derive the x
#' and y coordinates of the new points
#' @param together A boolean, if FALSE the eTraces plots will be generated for
#' each mapped point. Defaults to TRUE
#' @param new_colors The color(s) of the new points
#' @param jitter A boolean, if TRUE (default), it adds a random noise to place
#' the labels of new points a bit further apart. New runs lead to new positions
#'
#' @return The plot of the eTraces divided into two subplots, with the new points
#' highlighted
#' @export
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
#' edges_df <- data.frame('from' = edges_from, 'to' = edges_to, 'weight' = edges_weight,
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
#' mapped_obj <- mapNetwork(net, new_profiles)
#' map_eTrace(net, mapped_obj)
map_eTrace <- function(net,
                       mapped_obj,
                       genes = NULL,
                       score = NULL,
                       nRand=100,
                       main=NULL,
                       n_nearest = 15,
                       expression_enrichment = FALSE,
                       together = TRUE,
                       ylab = "score",
                       new_colors = NULL,
                       upper_colors = NULL,
                       lower_colors = NULL,
                       jitter = TRUE) {

  if(is.null(upper_colors)) {
    upper_colors <- net$Stages_color
  } else {
    if(length(upper_colors) != nrow(mapped_obj$new_cor)) {
      upper_colors <- rep(upper_colors, nrow(mapped_obj$new_cor))
    }
  }

  if(is.null(lower_colors)) {
    lower_colors <- net$SubClass_color
  } else {
    if(length(lower_colors) != nrow(mapped_obj$new_cor)) {
      lower_colors <- rep(lower_colors, nrow(mapped_obj$new_cor))
    }
  }

  names(upper_colors) <- net$Stages
  names(lower_colors) <- net$SubClass

  if(is.null(new_colors)) {
    derived_stage_color <- annotateMapping(net,
                                           new_cor = mapped_obj$new_cor,
                                           color_attr = 'Stages',
                                           col_vector = upper_colors,
                                           n_nearest = n_nearest)
    derived_stage_color <- upper_colors[derived_stage_color$Best.Annotation]

    derived_subclass_color <- annotateMapping(net,
                                              new_cor = mapped_obj$new_cor,
                                              color_attr = 'SubClass',
                                              col_vector = 'SubClass_color',
                                              n_nearest = n_nearest)
    derived_subclass_color <- lower_colors[derived_subclass_color$Best.Annotation]
  } else {
    if(length(new_colors) != ncol(mapped_obj$new_cor)) {
      new_colors <- rep(new_colors, ncol(mapped_obj$new_cor))
    }
    derived_stage_color <- new_colors
    derived_subclass_color <- new_colors
  }


  if(is.null(genes)) {
    genes <- rownames(net)
  }

  if(expression_enrichment) {
    eTrace <- get_eTrace(net = net, genes = genes, nRand = nRand)
    ylab <- "expression enrichment (z)"
  } else {
    if(is.null(score)) {
      if(!all(genes %in% rownames(net))) {
        message(paste0(genes[which(!genes %in% rownames(net))], collapse = '-'), ' not in `net`')
        genes <- genes[which(genes %in% rownames(net))]
        if(length(genes) == 0) {
          return('No gene in `net`')
        }
      }
      if(length(genes) == 1) {
        eTrace <- S4Vectors::List(z = SingleCellExperiment::logcounts(net)[genes,])
      } else {
        eTrace <- S4Vectors::List(z = Matrix::colMeans(SingleCellExperiment::logcounts(net)[genes,]))
      }
      ylab <- 'expression'
    } else {
      if(is.matrix(score)) {
        if(!all(genes %in% rownames(score))) {
          message(paste0(genes[which(!genes %in% rownames(score))], collapse = '-'), ' not in `score`')
          genes <- genes[which(genes %in% rownames(score))]
          if(length(genes) == 0) {
            return('No gene in `score`')
          }
        }
        if(length(genes) == 1) {
          eTrace <- S4Vectors::List(z = score[genes,])
        } else {
          eTrace <- S4Vectors::List(z = Matrix::colMeans(score[genes,]))
        }
      } else {
        eTrace <- S4Vectors::List(z = score)
      }
    }
  }

  xs <- seq(1,ncol(net))
  names(xs) <- colnames(net)

  ys <- eTrace$z
  names(ys) <- colnames(net)
  if(all(ys >= 0)) {
    col_zero_line <- NA
  } else {
    col_zero_line <- 'black'
  }

  derived_x <- apply(mapped_obj$new_cor, 2, function(i) {
    cors <- sort(i, decreasing = TRUE)[seq(1,n_nearest)]
    new_x <- stats::weighted.mean(xs[names(cors)], cors)
  })

  derived_y <- apply(mapped_obj$new_cor, 2, function(i) {
    cors <- sort(i, decreasing = TRUE)[seq(1,n_nearest)]
    new_y <- stats::weighted.mean(ys[names(cors)], cors)
  })

  sizes <- rep(1, length(ys))

  derived_sizes <- apply(mapped_obj$new_cor, 2, function(i) {
    new_size <- (mean(sort(i, decreasing = TRUE)[seq(1,n_nearest)])+1)
  })

  final_xs <- c(xs, derived_x)
  final_ys <- c(ys, derived_y)
  final_sizes <- c(sizes, derived_sizes)

  alphas <- c(rep(0.5, length(ys)), 1)

  nat_idx <- which(net$Stages == '6-early_infancy')[1]
  y_idx <- min(eTrace$z)+abs(min(eTrace$z))*0.05

  if(together) {
    graphics::par(mfrow = c(2, 1))

    graphics::par(mar = c(0, 5, 2, 2))
    plot(final_xs, final_ys,
         pch = c(rep(21, length(ys)), rep(22, length(derived_y))),
         bg  = scales::alpha(c(upper_colors, derived_stage_color),
                             c(rep(0.5, length(ys)), rep(1, length(derived_x)))),
         col = scales::alpha(c(upper_colors, rep("black", length(derived_x))),
                             c(rep(0.5, length(ys)), rep(1, length(derived_x)))),
         main = main, ylab = ylab,
         cex = final_sizes, xaxt = "n",
         lwd = c(rep(1, length(ys)), rep(2, length(derived_y))))

    graphics::abline(h = 0, lty = 2, lwd = 2, col = col_zero_line)
    graphics::abline(v = nat_idx, col = "darkgrey", lwd = 2, lty = 2)

    graphics::text(x = nat_idx + 5, y = y_idx, labels = "postnatal", pos = 4)
    graphics::text(x = nat_idx - 5, y = y_idx, labels = "prenatal", pos = 2)

    graphics::lines(stats::smooth.spline(seq(1, ncol(net)), eTrace$z, spar = 1),
                    col = "red", lwd = 2.5)

    usr <- graphics::par("usr")
    usr <- usr - usr/50

    if(jitter) {
      offset_x <- jitter(rep(0, length(derived_x)), factor = 200)
      offset_x[which(derived_x > stats::quantile(final_xs, 0.9) & offset_x > 0)] <- offset_x[which(derived_x > stats::quantile(final_xs, 0.9) & offset_x > 0)] * -1
      offset_x[which(derived_x < stats::quantile(final_xs, 0.1) & offset_x < 0)] <- offset_x[which(derived_x < stats::quantile(final_xs, 0.1) & offset_x < 0)] * -1
      offset_y <- jitter(rep(0, length(derived_y)), factor = 200)
      offset_y[which(derived_y > stats::quantile(final_ys, 0.9) & offset_y > 0)] <- offset_y[which(derived_y > stats::quantile(final_ys, 0.9) & offset_y > 0)] * -1
      offset_y[which(derived_y < stats::quantile(final_ys, 0.1) & offset_y < 0)] <- offset_y[which(derived_y < stats::quantile(final_ys, 0.1) & offset_y < 0)] * -1
    } else {
      offset_x <- 0
      offset_y <- 0
    }

    new_x <- pmin(pmax(derived_x + offset_x, usr[1]), usr[2])
    new_y <- pmin(pmax(derived_y + offset_y, usr[3]), usr[4])

    graphics::segments(x0 = derived_x, y0 = derived_y,
             x1 = new_x, y1 = new_y,
             col = "grey40", lty = 1)

    graphics::text(new_x, new_y, labels = colnames(mapped_obj$new_cor), cex = 0.8)

    graphics::par(mar = c(2.5, 5, 0.5, 2))
    plot(final_xs, final_ys,
         pch = c(rep(21, length(ys)), rep(22, length(derived_y))),
         bg  = scales::alpha(c(lower_colors, derived_subclass_color),
                             c(rep(0.5, length(ys)), rep(1, length(derived_x)))),
         col = scales::alpha(c(lower_colors, rep("black", length(derived_x))),
                             c(rep(0.5, length(ys)), rep(1, length(derived_x)))),
         main = "", ylab = ylab,
         cex = final_sizes,
         lwd = c(rep(1, length(ys)), rep(2, length(derived_y))))

    graphics::abline(h = 0, lty = 2, lwd = 2, col = col_zero_line)
    graphics::abline(v = nat_idx, col = "darkgrey", lwd = 2, lty = 2)

    graphics::text(x = nat_idx + 5, y = y_idx, labels = "postnatal", pos = 4)
    graphics::text(x = nat_idx - 5, y = y_idx, labels = "prenatal", pos = 2)

    graphics::lines(stats::smooth.spline(seq(1, ncol(net)), eTrace$z, spar = 1),
                    col = "red", lwd = 2.5)

    new_x2 <- pmin(pmax(derived_x + offset_x, usr[1]), usr[2])
    new_y2 <- pmin(pmax(derived_y + offset_y, usr[3]), usr[4])

    graphics::segments(x0 = derived_x, y0 = derived_y,
             x1 = new_x2, y1 = new_y2,
             col = "grey40", lty = 1)

    graphics::text(new_x2, new_y2, labels = colnames(mapped_obj$new_cor), cex = 0.8)
  } else {
    plots <- lapply(seq(1, length(derived_x)), function(i) {
      graphics::par(mfrow=c(2,1))
      graphics::par(mar=c(0,5,2,2))
      plot(c(xs, derived_x[i]), c(ys, derived_y[i]), pch=c(rep(21, length(ys)), 22),
           bg=scales::alpha(c(upper_colors, derived_stage_color), 0.5),
           col=scales::alpha(c(upper_colors, 'black'), c(rep(0.5, length(ys)), 1)),
           main=main,
           ylab=ylab,
           cex = final_sizes,
           xaxt = 'n',
           lwd = c(rep(1, length(ys)), 2))
      graphics::abline(h=0, lty=2, lwd=2, col = col_zero_line)
      graphics::abline(v=nat_idx, col = 'darkgrey', lwd = 2, lty = 2)
      graphics::text(x = nat_idx+5, y = y_idx, labels = 'postnatal', pos = 4)
      graphics::text(x = nat_idx-5, y = y_idx, labels = 'prenatal', pos = 2)
      graphics::lines(stats::smooth.spline(seq(1,ncol(net)),eTrace$z, spar = 1), col='red', lwd=2.5)
      graphics::text(x = derived_x[i], y = derived_y[i], labels = colnames(mapped_obj$new_cor)[i], pos = 3, cex = 0.8)

      graphics::par(mar=c(2.5,5,0.5,2))
      plot(c(xs, derived_x[i]), c(ys, derived_y[i]),
           pch=c(rep(21, length(ys)), 22),
           bg=scales::alpha(c(lower_colors, derived_subclass_color), 0.5),
           col=scales::alpha(c(lower_colors, 'black'), c(rep(0.5, length(ys)), 1)),
           main='',
           ylab=ylab,
           cex = final_sizes,
           xaxt = 'n',
           lwd = c(rep(1, length(ys)), 2))
      graphics::abline(h=0, lty=2, lwd=2, col = col_zero_line)
      graphics::abline(v=nat_idx, col = 'darkgrey', lwd = 2, lty = 2)
      graphics::text(x = nat_idx+5, y = y_idx, labels = 'postnatal', pos = 4)
      graphics::text(x = nat_idx-5, y = y_idx, labels = 'prenatal', pos = 2)
      graphics::lines(stats::smooth.spline(seq(1,ncol(net)),eTrace$z, spar = 1), col='red', lwd=2.5)
      graphics::text(x = derived_x[i], y = derived_y[i], labels = colnames(mapped_obj$new_cor)[i], pos = 3, cex = 0.8)
      p <- grDevices::recordPlot()
    })

    names(plots) <- colnames(mapped_obj$new_cor)

    return(plots)
  }
}
