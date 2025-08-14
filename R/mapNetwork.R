#' Map network
#'
#' @param net The reference network
#' @param new_profiles The new profiles to map
#' @param color_attr The color attribute in `net`, defaults to `SubClass`
#' @param label_attr The label attribute in `net`, defaults to `SubClass`
#' @param col_vector The color vector in `net`, defaults to `SubClass_color`
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
#' edges_df <- data.frame('from' = edges_from, 'to' = edges_to,
#' 'weight' = edges_weight,
#' 'from.x' = edges_from_x,
#' 'from.y' = edges_from_y,
#' 'to.x' = edges_to_x,
#' 'to.y' = edges_to_y)
#' net@metadata$network$edges <- edges_df
#' SummarizedExperiment::rowData(net)$informative <- sample(c(TRUE, FALSE),
#' size = nrow(net), replace = TRUE)
#' new_profiles <- matrix(sample(seq(1,10, length.out=10000), nrow(net)*10, replace = TRUE), ncol = 10)
#' rownames(new_profiles) <- rownames(net)
#' colnames(new_profiles) <- paste0('NewCol-', seq(1,10))
#' mapped_obj <- mapNetwork(net, new_profiles)
mapNetwork <- function(net,
                       new_profiles,
                       color_attr = 'SubClass',
                       label_attr = 'SubClass',
                       col_vector = NULL,
                       new_name = NULL,
                       n_nearest = 15, ...) {

  if(is.null(col_vector)) {
    col_vector <- paste0(color_attr, '_color')
  }

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
  names(col_vector) <- color_attr

  if(any(colnames(new_profiles) %in% color_attr)) {
    colnames(new_profiles)[which(colnames(new_profiles) %in% color_attr)] <- paste0('New-', colnames(new_profiles)[which(colnames(new_profiles) %in% color_attr)])
  }

  common_genes <- intersect(rownames(net)[SingleCellExperiment::rowData(net)$informative], rownames(new_profiles))
  new_cor <- stats::cor(t(apply(as.matrix(SingleCellExperiment::logcounts(net)[common_genes,]),
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
