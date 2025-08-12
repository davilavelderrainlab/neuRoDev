#' Map network
#'
#' @param net The reference network
#' @param new_profiles The new profiles to map
#' @param color_attr The color attribute in `net`, defaults to `SubClass`
#' @param label_attr The label attribute in `net`, defaults to `SubClass`
#' @param col_vector The color vector in `net`, defaults to `SubClass_colors`
#' @param new_name A vector to define the names of the new profiles. Defaults to
#' their colnames
#' @param n_nearest The number of nearest neighbors to consider for the
#' layout definition, mapping quality assessment and annotation. Defaults to 15
#' @param ... Other parameters for the `plotSameLayout` function
#'
#' @return An `S4Vectors` list containing the correlation matrix between the
#' profiles and the reference pseudobulks (`new_cor`), a data.frame containing
#' the edges information to connect the new profiles with the pseudobulks in the
#' network (`new_edges`), a matrix containing the layout information, including
#' the x and y coordinates (`new_layout`), the new plot (`new_plot`), a mapping
#' quality assessment (if mapping_quality = TRUE) (`mapping_quality`), and
#' the annotation scores (if annotate = TRUE) (`annotation`)
#' @export
#'
#' @examples
mapNetwork <- function(net,
                       new_profiles,
                       color_attr = 'SubClass',
                       label_attr = 'SubClass',
                       col_vector = 'SubClass_colors',
                       new_name = NULL,
                       n_nearest = 15, ...) {

  if(is.character(label_attr) && length(label_attr) == 1 &&
     label_attr %in% colnames(SingleCellExperiment::colData(net))) {
    label_attr <- SingleCellExperiment::colData(net)[,label_attr]
  }
  if(is.character(color_attr) && length(color_attr) == 1 &&
     color_attr %in% colnames(SingleCellExperiment::colData(net))) {
    color_attr <- SingleCellExperiment::colData(net)[,color_attr]
  }
  if(is.character(col_vector) && length(col_vector) == 1 &&
     col_vector %in% colnames(SingleCellExperiment::colData(net))) {
    col_vector <- SingleCellExperiment::colData(net)[,col_vector]
  }
  label_attr <- as.vector(label_attr)
  color_attr <- as.vector(color_attr)

  if(any(colnames(new_profiles) %in% color_attr)) {
    colnames(new_profiles)[which(colnames(new_profiles) %in% color_attr)] <- paste0('New-', colnames(new_profiles)[which(colnames(new_profiles) %in% color_attr)])
  }

  common_genes <- intersect(rownames(net)[SingleCellExperiment::rowData(net)$informative], rownames(new_profiles))
  new_cor <- stats::cor(t(apply(as.matrix(logcounts(net)[common_genes,]),
                                1, function(v) {(v-mean(v))/stats::sd(v)})),
                        new_profiles[common_genes,])

  out <- plotSameLayout(net,
                        new_cor = new_cor,
                        color_attr = color_attr,
                        label_attr = label_attr,
                        col_vector = col_vector,
                        new_name = new_name,
                        n_nearest = n_nearest, ...)

  return(out)

}
