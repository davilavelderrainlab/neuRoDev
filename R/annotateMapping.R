#' Add annotations to the mapped profiles
#'
#' @inheritParams plotSameLayout
#' @param compute_means A boolean, if TRUE the means of the neighbors scores are
#' considered instead of the sum (default)
#'
#' @return An `S4Vectors` list containing the global confidence of the mapping
#' annotations (`Global.Confidence`), the mapping confidence per mapped point
#' (`Mapping.Confidence`), the single annotations confidence per mapped point
#' (`Annotations.Confidence`), the best annotation per mapped point
#' (`Best.Annotation`), the highest confidence value per mapped point
#' (`Best.Confidence`), the annotations as a percentage (`Annotations.Perc`),
#' and a summary barplot (`Barplot`)
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
#' new_profiles <- matrix(sample(seq(1,10, length.out=10000), nrow(net)*10, replace = TRUE), ncol = 10)
#' rownames(new_profiles) <- rownames(net)
#' colnames(new_profiles) <- paste0('NewCol-', seq(1,10))
#' common_genes <- intersect(rownames(net)[SingleCellExperiment::rowData(net)$informative],
#' rownames(new_profiles))
#' new_cor <- stats::cor(t(apply(as.matrix(SingleCellExperiment::logcounts(net)[common_genes,]),
#' 1, function(v) {(v-mean(v))/stats::sd(v)})),new_profiles[common_genes,])
#' annotation <- annotateMapping(net, new_cor)
annotateMapping <- function(net,
                            new_cor,
                            color_attr = 'SubClass',
                            col_vector = NULL,
                            n_nearest = 15,
                            compute_means = FALSE,
                            order_names = FALSE) {

  if(is.null(col_vector)) {
    col_vector <- paste0(color_attr, '_color')
  }

  if(is.character(color_attr) && length(color_attr) == 1 && color_attr %in% colnames(SingleCellExperiment::colData(net))) {
    color_attr <- SingleCellExperiment::colData(net)[,color_attr]
  }
  if(is.character(col_vector) && length(col_vector) == 1 && col_vector %in% colnames(SingleCellExperiment::colData(net))) {
    col_vector <- SingleCellExperiment::colData(net)[,col_vector]
  }

  color_attr <- as.vector(color_attr)
  names(col_vector) <- color_attr

  idx_annotation <- which(apply(as.data.frame(SingleCellExperiment::colData(net)), 2, function(i) {
    all(as.vector(i) == color_attr)
  }))

  annotation_legend <- colnames(SingleCellExperiment::colData(net))[idx_annotation]

  n_nearest <- min(n_nearest, nrow(net)-1)

  annotation_tables <- list()

  for (c in colnames(new_cor)) {
    cors <- new_cor[,c]
    cors_filt <- sort(cors, decreasing = TRUE)[seq(1,n_nearest)]
    names(cors_filt) <- SingleCellExperiment::colData(net)[order(cors, decreasing = TRUE)[seq(1,n_nearest)],idx_annotation]
    annotations <- value_table(cors_filt,
                               perc = TRUE,
                               reciprocal = TRUE,
                               compute_means = compute_means)
    annotation_tables[[c]] <- annotations
  }

  mean_confidence <- apply(new_cor, 2,
                                function(i) {mean(sort(i, decreasing = TRUE)[seq(1,n_nearest)])})

  max_annotation <- unlist(lapply(annotation_tables, function(i) {
    names(i)[which.max(i)]
  }))
  max_value <- unlist(lapply(annotation_tables, function(i) {
    i[which.max(i)]
  }))
  names(max_value) <- max_annotation
  m <- get_matrix_from_list(annotation_tables)
  if(length(unique(as.vector(m))) == 1) {
    return("All clusters have the same scores, you should probably change what\n
           you map or the label attribute used.")
  }

  df <- reshape2::melt(m)
  df$Var1 <- as.character(df$Var1)
  if(order_names) {
    df$Var1 <- factor(df$Var1, levels = gtools::mixedsort(unique(df$Var1)))
  } else {
    df$Var1 <- factor(df$Var1, levels = unique(df$Var1))
  }
  df$Var2 <- factor(df$Var2, levels = gtools::mixedsort(unique(as.vector(df$Var2))))
  plot <- ggplot2::ggplot(data = df,
                          ggplot2::aes(x = Var1,
                                       y = value,
                                       fill = Var2)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(x = "",
                  y = "Values",
                  fill = annotation_legend,
                  title = '') +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_manual(values = col_vector[match(df$Var2, names(col_vector))]) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))

  confidence <- max_value*mean_confidence
  names(confidence) <- names(max_annotation)
  annotations_confidence <- lapply(seq(1,length(annotation_tables)), function(i) {
    annotation_tables[[i]]*mean_confidence[[i]]
  })
  names(annotations_confidence) <- names(annotation_tables)

  return(S4Vectors::SimpleList(Global.Confidence = mean(mean_confidence),
                               Mapping.Confidence = mean_confidence,
                               Annotations.Confidence = annotations_confidence,
                               Best.Annotation = max_annotation,
                               Best.Confidence = confidence,
                               Annotations.Perc = annotation_tables,
                               Barplot = plot))
}
