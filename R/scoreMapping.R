#' Mapping scores
#'
#' @inheritParams plotSameLayout
#'
#' @return An S4Vectors list containing the cluster accuracy `ClusterAccuracy`,
#' the label precision (`LabelPrecision`), the stage precision (`StagePrecision`),
#' the global accuracy (`GlobalAccuracy`), the global precision (`GlobalPrecision`),
#' the mapping score (`MappingScore`)
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
#' mapping_score <- scoreMapping(net, new_cor)
scoreMapping <- function(net,
                         new_cor,
                         label_attr = 'SubClass') {

  if(is.character(label_attr) && length(label_attr) == 1 && label_attr %in% colnames(SingleCellExperiment::colData(net))) {
    label_attr <- SingleCellExperiment::colData(net)[,label_attr]
  }
  label_attr <- as.vector(label_attr)

  mean_cors <- t(get_column_group_average(t(new_cor), interaction(label_attr, net$Stages, sep = 'and')))

  accuracy <- apply(mean_cors, 2, max)

  precision <- apply(mean_cors, 1, max)

  max_cor <- stats::cor(t(apply(as.matrix(SingleCellExperiment::logcounts(net)[rownames(net)[SingleCellExperiment::rowData(net)$informative],]), 1, function(v) {
    (v-mean(v))/stats::sd(v)
  })), as.matrix(SingleCellExperiment::logcounts(net)[rownames(net)[SingleCellExperiment::rowData(net)$informative],]))

  max_mean_cors <- t(get_column_group_average(t(max_cor), interaction(label_attr, net$Stages, sep = 'and')))

  max_accuracy <- mean(apply(max_mean_cors, 2, max))

  accuracy <- accuracy/max_accuracy

  max_precision <- apply(max_mean_cors, 1, max)

  precision <- precision/max_precision

  by_subclass_precision <- tapply(precision, unlist(lapply(strsplit(names(precision), 'and', fixed = TRUE), function(i) {i[1]})), mean)
  by_stage_precision <- tapply(precision, unlist(lapply(strsplit(names(precision), 'and', fixed = TRUE), function(i) {i[2]})), mean)

  global_accuracy <- mean(accuracy)

  accuracy[which(accuracy > 1)] <- 1

  global_precision <- mean(precision)

  mapping_score <- mean(c(global_accuracy, global_precision))

  return(S4Vectors::List(ClusterAccuracy = accuracy,
                         LabelPrecision = sort(by_subclass_precision),
                         StagePrecision = sort(by_stage_precision),
                         GlobalAccuracy = global_accuracy,
                         GlobalPrecision = global_precision,
                         MappingScore = mapping_score))

}
