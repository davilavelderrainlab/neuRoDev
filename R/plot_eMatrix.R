#' Title
#'
#' @inheritParams get_eMatrix
#' @param eMatrix If an eMatrix was already computed, it can be given here. If
#' NULL (default), it is computed.
#'
#' @return An Heatmap of the eMatrix
#' @export
#'
#' @examples
#' m <- matrix(sample(seq(1,10, length.out=10000), 15000*100, replace = TRUE), ncol = 100)
#' rownames(m) <- paste0('Gene-', seq(1,15000))
#' colnames(m) <- paste0('Col-', seq(1,100))
#' net <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = m))
#' net$SubClass <- rep(c('A', 'B', 'C', 'D'), each = 25)
#' net$Stages <- rep(c('S1', 'S2', 'S3', 'S4'), each = 25)
#' stages_palette <- c('S1' = 'pink', 'S2' = 'orange', 'S3' = 'violet', 'S4' = 'black')
#' net$Stages_color <- stages_palette[net$Stages]
#' plot_eMatrix(net = net, genes = rownames(net)[seq(1,5)])
plot_eMatrix <- function(net, genes, nRand=100, eMatrix = NULL) {

  if(is.null(eMatrix)) {
    eMatrix <- get_eMatrix(net = net, genes = genes, nRand = nRand)
  }

  z <- eMatrix$z

  col_fun <- circlize::colorRamp2(
    breaks = c(min(z, na.rm = TRUE), 0, max(z, na.rm = TRUE)),
    colors = c("blue", "white", "red")
  )

  stages_palette <- unique(net$Stages_color)
  names(stages_palette) <- unique(net$Stages)

  right_ha <- ComplexHeatmap::rowAnnotation(stage = names(stages_palette), col = list(stage = stages_palette), show_legend = FALSE, show_annotation_name = FALSE)
  ComplexHeatmap::Heatmap(z, cluster_rows = FALSE, cluster_columns = FALSE, col = col_fun, right_annotation = right_ha, show_row_names = FALSE, width = grid::unit(ncol(z)*4, 'mm'), height = grid::unit(nrow(z)*4, 'mm'), name = 'enrichment(z)')

}
