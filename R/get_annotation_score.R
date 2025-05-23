#' Get annotation score for mapping
#'
#' @param mapped_obj The result of mapData
#' @param reference_df The dataframe with the annotations per reference cluster
#' @param sub_idxs The indexes from which to select the annotations scores in
#' mapped_obj$NearestNeighborsAnnotation$Annotations. If NULL (default), they are
#' all kept.
#' @param annotation_attr Defaults to SubClass. The column name in reference_df in which
#' to find the cell annotation groups.
#'
#' @return A score per annotation considering all the mapped data points
#' @export
#'
#' @examples
#' set.seed(123)
#' S <- FC_signatures(matrix(runif(200,0,10), ncol = 10))
#' rownames(S) <- paste0('Gene-', seq(1, dim(S)[1]))
#' refS <- FC_signatures(matrix(runif(200,0.1,7), ncol = 10))
#' colnames(refS) <- paste0('Reference-', seq(1, dim(refS)[2]))
#' rownames(refS) <- paste0('Gene-', seq(1, dim(refS)[1]))
#' annotated_M <- reference_signatures_correlation(S, refS)
#' new_clusterS <- FC_signatures(matrix(runif(80,0,10), ncol = 4))
#' rownames(new_clusterS) <- paste0('Gene-', seq(1, dim(new_clusterS)[1]))
#' colnames(new_clusterS) <- paste0('New-', seq(1, dim(new_clusterS)[2]))
#' mapped_obj <- mapData(reference_df = annotated_M,
#' new_signatures = new_clusterS,
#' reference_signatures = refS,
#' color_attr = 'Best.Assignment')
#' get_annotation_score(mapped_obj, reference_df = annotated_M, annotation_attr = 'Best.Assignment')
get_annotation_score <- function(mapped_obj,
                                 reference_df,
                                 sub_idxs = NULL,
                                 annotation_attr = 'SubClass') {

  annotations <- mapped_obj$NearestNeighborsAnnotation$Annotations

  if(is.null(sub_idxs)) {
    sub_idxs <- seq(1, length(annotations))
  }

  annotations <- annotations[sub_idxs]

  annotation_score <- rowsum(as.matrix(unlist(annotations)),
                           group = unlist(lapply(strsplit(names(unlist(annotations)), '.', fixed = TRUE),
                                                 function(i) {paste0(i[seq(2, length(i))], collapse = '.')})))[names(table(unlist(lapply(strsplit(names(unlist(annotations)), '.', fixed = TRUE),
                                                                                                 function(i) {paste0(i[seq(2, length(i))], collapse = '.')})))),1]/table(unlist(lapply(strsplit(names(unlist(annotations)), '.', fixed = TRUE),
                                                                                                                                               function(i) {i[2]})))

  for(annotation in unique(reference_df[,annotation_attr])) {

    if(!annotation %in% names(annotation_score)) {
      previous_names <- names(annotation_score)
      annotation_score <- c(annotation_score, 0)
      names(annotation_score) <- c(previous_names, annotation)
    }

  }

  annotation_score <- annotation_score[unique(reference_df[,annotation_attr])]

  return(annotation_score)

}
