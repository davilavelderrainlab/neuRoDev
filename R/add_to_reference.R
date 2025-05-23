#' Add new clusters to an annotated reference
#'
#'  Given an annotated reference dataframe (reference_df),
#'  a correlation dataframe (signatures_cor), as defined with
#'  reference_signature_correlation, an annotation label (and possibly
#'  label_attr), it maps the new clusters on the annotated reference
#'  providing a UMAP object and a UMAP plot (multiple plots if focus = TRUE and
#'  label_attr is given, with each plot that zooms-in each group
#'  defined in color_attr, subdividing the clusters based on
#'  label_attr). The default method of clustering the UMAP is
#'  Louvain ('louvain'), but it could be also Walktrap ('walktrap') and
#'  Leiden ('leiden'). A title to the UMAP is given if title is provided.
#'  You can provide color palettes for the UMAP for both the annotation and the sub annotation
#' @inheritParams umap_plot_same_layout
#' @param color_attr The labels to use for the annotation/coloring of the plot
#' @param label_attr In case annotation labels different than the ones of the colors are wanted, this specifies the labels to show
#' @param no_label A boolean variable. If TRUE, no label is added
#' @param ... Additional parameters for `umap_plot_same_layout`
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
#' new_M <- reference_signatures_correlation(new_clusterS, refS)
#' add_to_reference(annotated_M,
#' new_M,
#' annotated_M$`Best.Assignment`)
add_to_reference <- function(reference_df,
                             signatures_cor,
                             color_attr,
                             label_attr=NULL,
                             new_name=NULL,
                             no_label=FALSE,
                             only_new_points=FALSE,
                             ...) {

  if(only_new_points && no_label) {
    message('Only the new points will be labeled in the plot,
                even though no_label was true. If no label is wanted,
                set only_new_points as FALSE.')
  }

  if(is.character(color_attr) && length(color_attr) == 1 && color_attr %in% colnames(reference_df)) {
    color_attr <- reference_df[,color_attr]
  }

  if(is.character(label_attr) && length(label_attr) == 1 && label_attr %in% colnames(reference_df)) {
    label_attr <- reference_df[,label_attr]
  }

  if(is.null(new_name)) {new_name <- rownames(signatures_cor)}

  if(any(new_name %in% c(color_attr, label_attr))) {

    new_name[which(new_name %in% c(color_attr, label_attr))] <- paste0('New-', new_name[which(new_name %in% c(color_attr, label_attr))])

  }

  if (any(rownames(signatures_cor) %in% c(color_attr, label_attr))) {

    rownames(signatures_cor)[which(rownames(signatures_cor) %in% c(color_attr, label_attr))] <- paste0("New-", rownames(signatures_cor)[which(rownames(signatures_cor) %in% c(color_attr, label_attr))])

  }

  if(length(new_name) == 1) {
    color_attr <- c(color_attr, rep(new_name, dim(signatures_cor)[1]))
  } else {
    color_attr <- c(color_attr, new_name)
  }

  if(no_label & !is.null(label_attr)) {
    message('no_label is set to TRUE but a label_attr is given. No labels will be shown. To avoid this, set no_label to FALSE')
    label_attr <- NULL
  } else if(!no_label & is.null(label_attr)) {
    label_attr <- color_attr
  }

  action_umap <- umap_plot_same_layout(reference_df = reference_df,
                                       signatures_cor = signatures_cor,
                                       color_attr = color_attr,
                                       label_attr = label_attr,
                                       new_name = new_name,
                                       only_new_points = only_new_points, ...)

  return(action_umap)

}
