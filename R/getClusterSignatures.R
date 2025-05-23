#' Signatures Initialization
#'
#' @param expression_matrix A single-cell RNAseq matrix or a list of matrices
#' @param split A vector that defines how to split a expression_matrix into multiple ones
#' @param subsample A boolean variable that defines if each of the provided
#' expression_matrixs has to be proportionally subsampled to the same number of cells
#' @param PCs_to_use The number of Principal Components to consider for
#' FindNeighbors of Seurat
#' @param resolution The resolution to consider for FindClusters of Seurat
#' @param clusters_mv A vector that substitutes the clustering
#' (ex. celltypes if already annotated)
#' @param Ntotal The total number of cells to consider in the subsampling
#' @param BOGs A boolean variable to define if map_clusters_to_representation has
#' to compute also the BOGs
#' @param verbose A boolean variable to define if the steps of Seurat will be
#' printed or not
#'
#' @return The result of map_clusters_to_representation run on each expression_matrix after
#' splitting, subsampling, clustering if/when required
#' @export
#'
#' @examples
#' set.seed(123)
#' M <- matrix(runif(2000000,0,100), ncol=20000)
#' rownames(M) <- paste0('Gene-', seq(1, dim(M)[1]))
#' colnames(M) <- paste0('Cell-', seq(1, dim(M)[2]))
#' getClusterSignatures(M, resolution=2)
getClusterSignatures <- function(expression_matrix,
                                 split=NULL,
                                 subsample=FALSE,
                                 PCs_to_use=NULL,
                                 resolution=1,
                                 clusters_mv=NULL,
                                 Ntotal=2000,
                                 BOGs=FALSE,
                                 verbose=FALSE) {

  if(!is.list(expression_matrix)) {
    expression_matrix <- Matrix::Matrix(expression_matrix, sparse = TRUE)
    expression_matrix <- list(expression_matrix)
    if(!is.null(split)) {
      split <- list(split)
    }
    if(!is.null(clusters_mv)) {
      clusters_mv <- list(clusters_mv)
    }
  }

  if(!is.null(split)) {
    expression_matrix <- lapply(seq_len(length(expression_matrix)), function(i) {
      matrix_split(expression_matrix[[i]], split[[i]])
    })
    expression_matrix <- unlist(expression_matrix, recursive = FALSE)
  }

  expression_matrix <- lapply(expression_matrix, function(i) {Matrix::Matrix(i, sparse = TRUE)})

  seurat_objects <- lapply(expression_matrix,
                           run_seurat_clustering,
                           PCs_to_use = PCs_to_use,
                           resolution = resolution,
                           verbose = verbose)

  if(is.null(clusters_mv)) {

    mvs <- lapply(seurat_objects, function(s) {s$seurat_clusters})

    if(subsample) {

      Ntotal <- min(unlist(lapply(expression_matrix, function(i) {dim(i)[2]})), Ntotal)

      new_expression_matrix <- lapply(seq_len(length(expression_matrix)), function(i) {
        Proportional_sampling(expression_matrix[[i]], mvs[[i]], Ntotal=Ntotal)
      })

      new_seurat_objects <- lapply(new_expression_matrix,
                                   run_seurat_clustering,
                                   PCs_to_use = PCs_to_use,
                                   resolution = resolution,
                                   verbose = verbose)

      new_mvs <- lapply(new_seurat_objects, function(s) {s$seurat_clusters})

      summary <- lapply(seq_len(length(expression_matrix)), function(i) {
        if(length(unique(mvs[[i]])) == 1) {
          return(S4Vectors::List(group_counts = mvs[[i]],
                                 P = as.matrix(Matrix::rowMeans(new_seurat_objects[[i]]@assays$RNA$data)),
                                 S = as.matrix(Matrix::rowMeans(new_seurat_objects[[i]]@assays$RNA$data)),
                                 BOG = NULL,
                                 DE_out = NULL))
        }
        map_clusters_to_representation(new_seurat_objects[[i]]@assays$RNA$data,
                                       new_mvs[[i]],
                                       BOGs = BOGs)
      })

      summary <- lapply(seq(1,length(summary)), function(i) {
        s <- summary[[i]]
        s$Membership <- new_mvs[[i]]
        return(s)
      })

    } else {

      summary <- lapply(seq_len(length(expression_matrix)), function(i) {
        if(length(unique(mvs[[i]])) == 1) {
          return(S4Vectors::List(group_counts = mvs[[i]],
                                 P = as.matrix(Matrix::rowMeans(seurat_objects[[i]]@assays$RNA$data)),
                                 S = as.matrix(Matrix::rowMeans(seurat_objects[[i]]@assays$RNA$data)),
                                 BOG = NULL,
                                 DE_out = NULL))
        }
        map_clusters_to_representation(seurat_objects[[i]]@assays$RNA$data,
                                       mvs[[i]],
                                       BOGs = BOGs)
      })

      summary <- lapply(seq(1,length(summary)), function(i) {
        s <- summary[[i]]
        s$Membership <- mvs[[i]]
        return(s)
      })

      new_seurat_objects <- NULL
      new_mvs <- NULL

    }

  } else {

    if(subsample) {

      Ntotal <- min(unlist(lapply(expression_matrix, function(i) {dim(i)[2]})), Ntotal)

      new_expression_matrix <- lapply(seq_len(length(expression_matrix)), function(i) {
        Proportional_sampling(expression_matrix[[i]], clusters_mv[[i]], Ntotal=Ntotal)
      })

      new_seurat_objects <- lapply(new_expression_matrix,
                                   run_seurat_clustering,
                                   PCs_to_use = PCs_to_use,
                                   resolution = resolution,
                                   verbose = verbose)

      new_clusters_mv <- lapply(seq_len(length(expression_matrix)), function(i) {
        idxs <- Proportional_sampling(expression_matrix[[i]],
                                      clusters_mv[[i]],
                                      Ntotal=Ntotal,
                                      ReturnSCE = FALSE)$idx
        clusters_mv[[i]][idxs]
      })

      summary <- lapply(seq_len(length(new_expression_matrix)), function(i) {
        if(length(unique(mvs[[i]])) == 1) {
          return(S4Vectors::List(group_counts = mvs[[i]],
                                 P = as.matrix(Matrix::rowMeans(new_seurat_objects[[i]]@assays$RNA$data)),
                                 S = as.matrix(Matrix::rowMeans(new_seurat_objects[[i]]@assays$RNA$data)),
                                 BOG = NULL,
                                 DE_out = NULL))
        }

        map_clusters_to_representation(new_seurat_objects[[i]]@assays$RNA$data,
                                       new_clusters_mv[[i]],
                                       BOGs = BOGs)
      })

      summary <- lapply(seq(1,length(summary)), function(i) {
        s <- summary[[i]]
        s$Membership <- new_clusters_mv[[i]]
        return(s)
      })

    } else {

      summary <- lapply(seq_len(length(expression_matrix)), function(i) {
        if(length(unique(clusters_mv[[i]])) == 1) {
          return(S4Vectors::List(group_counts = clusters_mv[[i]],
                                 P = as.matrix(Matrix::rowMeans(seurat_objects[[i]]@assays$RNA$data)),
                                 S = as.matrix(Matrix::rowMeans(seurat_objects[[i]]@assays$RNA$data)),
                                 BOG = NULL,
                                 DE_out = NULL))
        }
        map_clusters_to_representation(seurat_objects[[i]]@assays$RNA$data,
                                       clusters_mv[[i]],
                                       BOGs = BOGs)
      })

      summary <- lapply(seq(1,length(summary)), function(i) {
        s <- summary[[i]]
        s$Membership <- clusters_mv[[i]]
        return(s)
      })

      new_seurat_objects <- NULL

    }
  }

  if(is.null(new_seurat_objects)) {
    new_seurat_objects <- seurat_objects
  }

  if(length(summary) == 1) {
    summary <- summary[[1]]
  }

  if(length(new_seurat_objects) == 1) {
    new_seurat_objects <- new_seurat_objects[[1]]
  }

  return(S4Vectors::List('SeuratObject' = new_seurat_objects, 'Summary' = summary))

}
