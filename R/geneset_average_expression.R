#' Get average expression for genesets or geneset groups
#'
#' @param geneset_df The DF that specifies the wanted genesets
#' @param pseudobulk The wanted pseudobulk expression
#' @param genesets The column in which the genesets are specified in
#' `geneset_df`
#' @param celltypes The column in which the celltypes are specified in
#' `geneset_df`
#' @param genes The column in which the wanted genes are specified in
#' `geneset_df`
#'
#' @return The average expression of each pseudobulk column for each geneset
#' @export
#'
#' @examples
#' anno_df <- S4Vectors::DataFrame('AnnotationGroup' = paste0('AnnoGroup-',
#' seq(1, 10)),
#' 'Class' = rep('Class1', 10),
#' 'TopGenesByGroup' = rep('Gene.1-Gene.2', 10))
#' pseudo <- SingleCellExperiment::SingleCellExperiment(
#' assays = list(logcounts = matrix(sample(seq(1,100),
#' 4*10,
#' replace = TRUE),
#' nrow = 4)))
#' rownames(pseudo) <- paste0('Gene.', seq(1, nrow(pseudo)))
#' colnames(pseudo) <- paste0('Column-', seq(1, ncol(pseudo)))
#' geneset_average_expression(anno_df, pseudo)
geneset_average_expression <- function(geneset_df,
                                       pseudobulk,
                                       genesets = 'AnnotationGroup',
                                       celltypes = 'Class',
                                       genes = 'TopGenesByGroup') {

  if(is(pseudobulk, 'SingleCellExperiment')) {
    pseudobulk <- SingleCellExperiment::logcounts(pseudobulk)
  }

  geneset_matrix <- do.call(rbind, lapply(unique(interaction(geneset_df[,celltypes],
                                                                geneset_df[,genesets])), function(i) {

                                                                  i <- as.vector(i)
                                                                  ct <- unlist(lapply(strsplit(i, '.', fixed = TRUE), function(i) {i[1]}))
                                                                  ag <- unlist(lapply(strsplit(i, '.', fixed = TRUE), function(i) {i[2]}))

                                                                  sub_anno_df <- geneset_df[which(geneset_df[,genesets] == ag & geneset_df[,celltypes] == ct),]
                                                                  all_genes <- unlist(strsplit(sub_anno_df[,genes], '-', fixed = TRUE))

                                                                  sub_exp <- SingleCellExperiment::logcounts(pseudobulk)[all_genes,]

                                                                  mean_exp <- Matrix::colMeans(sub_exp)

                                                                }))

  return(geneset_matrix)

}
