#' Get the expression enrichment Trace
#'
#' @param net The reference network
#' @param genes Genes of which to compute the expression enrichment
#' @param nRand The number of random sampling to use. If NULL, a nRand
#' proportional to the number of genes is computed (the higher the number of
#' genes, the lower the number of nRand computed, to a minimum of 200).
#'
#' @return An S4Vectors List with the averaged scaled logcounts of the genes (`obs`),
#' the z score of the enrichment (`z`), and the pvalue (`p`)
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
#' eTrace <- get_eTrace(net, 'Gene-1')
get_eTrace <- function(net,
                       genes,
                       nRand = 500L) {

  x <- SingleCellExperiment::logcounts(net)
  rn <- rownames(x)
  present <- genes %in% rn

  if(!all(present)) {
    message(paste0(genes[!present], collapse = "-"), " not in `net`")
    genes <- genes[present]
    if(length(genes) == 0L) {
      return("No gene in `net`")
    }
  }

  genes <- unique(genes)

  nc <- ncol(x)
  if(nc < 2L) {
    stop("Need at least 2 columns in `logcounts(net)` to compute row SDs.")
  }

  rm <- Matrix::rowMeans(x)
  xc <- sweep(x, 1L, rm, "-", check.margin = FALSE)
  rsd <- sqrt(Matrix::rowSums(xc^2) / (nc - 1L))

  ok_row <- is.finite(rsd) & (rsd > 0)
  xc <- xc[ok_row, , drop = FALSE]
  rsd <- rsd[ok_row]

  scaled_exp <- sweep(xc, 1L, rsd, "/", check.margin = FALSE)

  genes <- intersect(genes, rownames(scaled_exp))
  if(length(genes) == 0L) {
    return("No non-constant gene in `net` after scaling")
  }

  y <- Matrix::colMeans(scaled_exp[genes, , drop = FALSE])

  keep <- !(rownames(scaled_exp) %in% genes)
  bg <- scaled_exp[keep, , drop = FALSE]

  n_bg <- nrow(bg)
  leng <- length(genes)

  if(is.null(nRand)) {
    fnr <- nrow(scaled_exp)
    nRand <- floor(fnr/((round(fnr/200)*leng)/(leng+50)))
  }

  if(n_bg == 0L) {
    stop("No background rows left after removing `genes`.")
  }

  if(leng > 1L) {
    if(n_bg < leng) {
      stop("Not enough background rows to sample random groups of size length(genes).")
    }
    n_done <- 0L
    mean_acc <- numeric(nc)
    m2_acc <- numeric(nc)

    starts <- seq.int(1L, nRand, by = 100L)

    for(s in starts) {
      e <- min(s + 100L - 1L, nRand)
      m <- e - s + 1L

      idx_chunk <- replicate(m, sample.int(n_bg, leng, replace = FALSE))

      chunk_means <- vapply(
        seq_len(m),
        function(j) Matrix::colMeans(bg[idx_chunk[, j], , drop = FALSE]),
        numeric(nc)
      )

      for(j in seq_len(m)) {
        n_done <- n_done + 1L
        xj <- chunk_means[, j]
        delta <- xj - mean_acc
        mean_acc <- mean_acc + delta / n_done
        delta2 <- xj - mean_acc
        m2_acc <- m2_acc + delta * delta2
      }
    }

    cm <- mean_acc
    csd <- sqrt(m2_acc / (n_done - 1L))

  } else {
    if(n_bg < 2L) {
      stop("Need at least 2 background rows to compute column SDs.")
    }

    cm <- Matrix::colMeans(bg)
    bgc <- sweep(bg, 2L, cm, "-", check.margin = FALSE)
    csd <- sqrt(Matrix::colSums(bgc ^ 2) / (n_bg - 1L))
  }

  csd[!is.finite(csd) | csd == 0] <- NA_real_

  z <- (y - cm) / csd
  p <- 2 * stats::pnorm(-abs(z))

  S4Vectors::List(obs = y,
                  z = z,
                  p = p)

}
