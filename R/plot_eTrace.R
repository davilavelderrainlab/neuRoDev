#' Plot the eTraces
#'
#' @inheritParams get_eTrace
#' @param score A score to plot on the y axis
#' @param expression_enrichment A boolean, if TRUE the expression enrichment is
#' computed and plotted (`get_eTrace` function). Defaults to FALSE
#' @param main The title of the plot
#' @param upper_colors The colors for the upper plot
#' @param lower_colors The colors for the lower plot
#' @param ylab The ylab to add. Defaults to `score`
#'
#' @return The plot of the eTraces divided into two subplots
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
#' plot_eTrace(net, 'Gene-1')
plot_eTrace <- function(net,
                        genes = NULL,
                        score = NULL,
                        expression_enrichment = FALSE,
                        main=NULL,
                        upper_colors = NULL,
                        lower_colors = NULL,
                        nRand = 100,
                        ylab = "score") {

  if(is.null(upper_colors)) {
    upper_colors <- net$Stages_color
  }

  if(is.null(lower_colors)) {
    lower_colors <- net$SubClass_color
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

  col_zero_line <- ifelse(all(eTrace$z >= 0), NA, 'black')

  x <- seq(1,ncol(net))
  nat_idx <- which(net$Stages == '6-early_infancy')[1]
  y_idx <- min(eTrace$z)+abs(min(eTrace$z))*0.05

  graphics::par(mfrow=c(2,1))
  graphics::par(mar=c(0,5,2,2))
  plot(eTrace$z, pch=21, bg=upper_colors, main=main, ylab=ylab, xaxt = 'n')
  graphics::abline(h=0, lty=2, lwd=2, col = col_zero_line)
  graphics::abline(v=nat_idx, col = 'darkgrey', lwd = 2, lty = 2)
  graphics::text(x = nat_idx+5, y = y_idx, labels = 'postnatal', pos = 4)
  graphics::text(x = nat_idx-5, y = y_idx, labels = 'prenatal', pos = 2)
  graphics::lines(stats::smooth.spline(x,eTrace$z, spar = 1), col='red', lwd=2.5)
  graphics::par(mar=c(2.5,5,0.5,2))
  plot(eTrace$z, pch=21, bg=lower_colors, ylab=ylab)
  graphics::abline(h=0, lty=2, lwd=2, col = col_zero_line)
  graphics::abline(v=nat_idx, col = 'darkgrey', lwd = 2, lty = 2)
  graphics::text(x = nat_idx+5, y = y_idx, labels = 'postnatal', pos = 4)
  graphics::text(x = nat_idx-5, y = y_idx, labels = 'prenatal', pos = 2)
  graphics::lines(stats::smooth.spline(x,eTrace$z, spar = 1), col='red', lwd=2.5)
}
