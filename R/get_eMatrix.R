#' Get the expression enrichment matrix
#'
#' @param net The reference network
#' @param genes Genes of which to compute the expression enrichment
#' @param nRand The number of random sampling to use (defaults to 100)
#'
#' @return A list with z-scores and p-values matrices divided into
#' subclass-stage pairs
#' @export
#'
#' @examples
#' m <- matrix(sample(seq(1,10, length.out=10000), 15000*100, replace = TRUE), ncol = 100)
#' rownames(m) <- paste0('Gene-', seq(1,15000))
#' colnames(m) <- paste0('Col-', seq(1,100))
#' net <- SingleCellExperiment::SingleCellExperiment(assays = list(logcounts = m))
#' net$SubClass <- rep(c('A', 'B', 'C', 'D'), each = 25)
#' net$Stages <- rep(c('S1', 'S2', 'S3', 'S4'), each = 25)
#' emat <- get_eMatrix(net = net, genes = rownames(net)[seq(1,5)])
get_eMatrix <- function(net, genes, nRand = NULL) {

  if(is.null(nRand)) {
    nr <- nrow(net)
    lgn <- length(genes)
    nRand <- floor(nr/((round(nr/200)*lgn)/(lgn+50)))
  }

  nRand <- as.integer(nRand)
  if (nRand < 1L) stop("nRand must be >= 1.")

  X <- SummarizedExperiment::assay(net, "logcounts")
  gene_names <- rownames(net)

  test_genes <- intersect(genes, gene_names)
  if (!length(test_genes)) stop("No genes found in net rownames().")

  k_rand <- length(test_genes)

  if (k_rand > nrow(X)) stop("length(genes) cannot exceed nrow(net).")

  test_idx <- match(test_genes, gene_names)

  if(length(test_genes)==1) {
    nRand <- nrow(X)
    rand_idx <- lapply(seq(1,nRand), function(i) {i})
  } else {
    rand_idx <- replicate(
      nRand,
      sample.int(n = nrow(X), size = k_rand, replace = FALSE),
      simplify = FALSE
    )
  }

  all_sets <- c(list(test_idx), rand_idx)
  nSets <- length(all_sets)

  union_idx <- sort.int(unique(unlist(all_sets, use.names = FALSE)))

  W <- Matrix::sparseMatrix(
    i = match(unlist(all_sets, use.names = FALSE), union_idx),
    j = rep.int(seq_len(nSets), lengths(all_sets)),
    x = 1,
    dims = c(length(union_idx), nSets)
  )

  stages   <- SummarizedExperiment::colData(net)[["Stages"]]
  subclass <- SummarizedExperiment::colData(net)[["SubClass"]]

  stage_levels <- sort(unique(stages))
  sub_levels   <- sort(unique(subclass))

  nStage <- length(stage_levels)
  nSub   <- length(sub_levels)

  res <- array(
    0,
    dim = c(nSets, nStage, nSub),
    dimnames = list(NULL, as.character(stage_levels), as.character(sub_levels))
  )

  stage_idx <- lapply(stage_levels, function(st) which(stages == st))
  names(stage_idx) <- as.character(stage_levels)

  stage_sub_idx <- lapply(stage_idx, function(idx) {
    out <- lapply(sub_levels, function(sb) which(subclass[idx] == sb))
    names(out) <- as.character(sub_levels)
    out
  })

  for (s in seq_along(stage_levels)) {
    idx <- stage_idx[[s]]
    if (!length(idx)) next

    M <- as.matrix(X[union_idx, idx, drop = FALSE])

    mu  <- matrixStats::rowMeans2(M)
    sdv <- matrixStats::rowSds(M)
    sdv[!is.finite(sdv) | sdv == 0] <- NA_real_

    Z <- sweep(M, 1L, mu, "-", check.margin = FALSE)
    Z <- sweep(Z, 1L, sdv, "/", check.margin = FALSE)

    ok <- !is.na(Z)
    Z[!ok] <- 0
    ok_num <- ok * 1

    num <- Matrix::crossprod(W, Z)
    den <- Matrix::crossprod(W, ok_num)
    score <- as.matrix(num / den)

    out_stage <- matrix(
      0,
      nrow = nSets,
      ncol = nSub,
      dimnames = list(NULL, as.character(sub_levels))
    )

    for (k in seq_along(sub_levels)) {
      cols <- stage_sub_idx[[s]][[k]]
      if (!length(cols)) next

      out_stage[, k] <- matrixStats::rowMedians(
        score[, cols, drop = FALSE],
        na.rm = FALSE
      )
    }

    out_stage[is.na(out_stage)] <- 0
    res[, s, ] <- out_stage
  }

  Test <- matrix(
    res[1, , , drop = FALSE],
    nrow = nStage,
    ncol = nSub,
    dimnames = list(as.character(stage_levels), as.character(sub_levels))
  )

  rand_mat <- matrix(res[-1, , , drop = FALSE], nrow = nRand)

  Ave <- matrix(
    colMeans(rand_mat),
    nrow = nStage,
    ncol = nSub,
    dimnames = dimnames(Test)
  )

  SDs <- matrix(
    matrixStats::colSds(rand_mat),
    nrow = nStage,
    ncol = nSub,
    dimnames = dimnames(Test)
  )

  z <- (Test - Ave) / SDs
  p <- 2 * stats::pnorm(-abs(z))

  S4Vectors::List(z = z, p = p)
}
