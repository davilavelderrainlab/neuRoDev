#' Run standard Seurat pipeline given a matrix
#'
#' @param expression_matrix An expression matrix
#' @param PCs_to_use The number of Principal Components to consider for
#' FindNeighbors of Seurat
#' @param verbose A boolean variable to define if the steps of Seurat will be
#' printed or not
#' @param resolution The resolution to consider for FindClusters of Seurat
#'
#' @return A seurat object with clusters, UMAP, PCA and the normalized counts
#' @export
#'
#' @examples
#' set.seed(123)
#' M <- matrix(runif(200000,0,100), ncol=2000)
#' rownames(M) <- paste0('Gene-', seq(1, dim(M)[1]))
#' colnames(M) <- paste0('Cell-', seq(1, dim(M)[2]))
#' run_seurat_clustering(M)
run_seurat_clustering <- function(expression_matrix,
                                  PCs_to_use=NULL,
                                  verbose=FALSE,
                                  resolution=1) {

  if(is.null(colnames(expression_matrix))) {
    colnames(expression_matrix) <- paste0("c", seq(1, ncol(expression_matrix)))
  }

  s <- Seurat::CreateSeuratObject(counts = expression_matrix)
  s <- Seurat::NormalizeData(object = s,
                             verbose = verbose)

  span <- 0.3
  span_x <- 0.3
  warns <- TRUE

  while(warns) {
    tryCatch(Seurat::FindVariableFeatures(object = s,
                                          verbose = verbose,
                                          span = span), warning = function(x) {
                                            assign('span_x', span + 0.1, envir = '.GlobalEnv')
                                          })
    if(span_x != span) {
      span <- span_x
    } else {
      warns <- FALSE
    }
  }

  s <- Seurat::FindVariableFeatures(object = s,
                                    verbose = verbose,
                                    span = span)
  s <- Seurat::ScaleData(object = s,
                         features = Seurat::VariableFeatures(s),
                         verbose = verbose)
  numPCs <- min(50, dim(expression_matrix)[2]-1)
  s <- Seurat::RunPCA(object = s,
                      features = Seurat::VariableFeatures(s),
                      npcs = numPCs,
                      approx = FALSE,
                      verbose = verbose)

  if(is.null(PCs_to_use)) {
    PCs_to_use <- (s$pca@stdev)^2
    PCs_to_use <- PCs_to_use/sum(PCs_to_use)
    PCs_to_use <- cumsum(PCs_to_use)
    PCs_to_use <- min(which(PCs_to_use >= 0.75))
  }

  PCs_to_use <- min(numPCs, PCs_to_use)

  s <- Seurat::FindNeighbors(object = s,
                             dims = seq_len(PCs_to_use),
                             verbose = verbose)

  s <- Seurat::FindClusters(s, resolution = resolution, verbose = verbose)
  s$seurat_clusters <- paste0("c", s$seurat_clusters)
  nn <- min(30, dim(expression_matrix)[2] - 1)

  if(nn < 2) {
    return(s)
  }

  s <- Seurat::RunUMAP(s,
                       dims = seq(1,10),
                       verbose = verbose,
                       n.neighbors = nn,
                       umap.method = "uwot",
                       metric = "cosine")

  return(s)
}
