#' Map new data onto a reference network
#'
#' @param annotated_reference The annotated reference correlation data frame
#' @param new_profiles The profiles that have to be mapped (can directly provide
#' signatures)
#' @param new_signatures The signatures that have to be mapped (can be computed
#' from profiles)
#' @param reference_signatures The reference signatures used to build the
#' reference network
#' @param reference_layout The reference network layout object
#' @param color_attr The annotation label wanted to be displayed and used
#' for the annotation of the new points
#' @param label_attr The labels to be used in the UMAP. If NULL, the labels of
#' color_attr will be used, otherwise color_attr will be used to
#' color the clusters and label_attr to show the labels.
#' @param col_vector The color vector palette. Should be a named vector where
#' names correspond to the unique values in the color_attr vector.
#' @param ... Additional parameters given to `add_to_annotated_reference`
#'
#' @return The umap object and plot before and after adding the new clusters,
#' plus mapping quality metrics and nearest neighbor annotation
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
#' mapData(annotated_reference = annotated_M,
#' new_signatures = new_clusterS,
#' reference_signatures = refS,
#' color_attr = 'Best.Assignment')
mapData <- function(annotated_reference,
                    new_profiles=NULL,
                    new_signatures=NULL,
                    reference_signatures,
                    color_attr,
                    label_attr=NULL,
                    reference_layout=NULL,
                    col_vector=NULL,
                    ...) {

  if(is.null(new_signatures) & is.null(new_profiles)) {
    return('Error: either signatures or profiles need to be given.')
  }

  if(!is.null(new_signatures) & !is.null(new_profiles)) {
    message('Both signatures and profiles are given. Only the first will be used.')
  }

  if(is.null(new_signatures)) {
    new_signatures <- FC_signatures(new_profiles)
  }

  mean_sig <- mean(unlist(new_signatures) > 1)
  if(mean_sig < 0.05) {
    message('The signatures given have less then 5% (', round(mean_sig*100, digits = 1), '%)
    of the genes with significant variability in the samples (fold-change > 2).
    The mapping results should be interpreted with caution, as the small
    differences will likely be amplified. This may also happen if there are multiple
    clusters of similar subtypes.')
  }
  if(mean_sig > 0.20) {
    message('The signatures given have more then 20% (', round(mean_sig*100, digits = 1), '%)
    of the genes with significant variability in the samples (fold-change > 2).
    Check if the signatures were computed on normalized data for better results')
  }

  new_cor <- reference_signatures_correlation(S = new_signatures,
                                              refS = reference_signatures,
                                              warn = FALSE)

  out <- add_to_annotated_reference(annotated_reference = annotated_reference,
                                    signatures_cor = new_cor,
                                    color_attr = color_attr,
                                    label_attr = label_attr,
                                    new_name = rownames(new_cor),
                                    umap_obj = reference_layout,
                                    col_vector = col_vector, ...)

  return(out)

}
