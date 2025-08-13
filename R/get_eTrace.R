#' Get the expression enrichment Trace
#'
#' @param net The reference network
#' @param genes Genes of which to compute the expression enrichment
#' @param nRand The number of random sampling to use (defaults to 100)
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
                       nRand=100) {

  if(!all(genes %in% rownames(net))) {
    message(genes[which(!genes %in% rownames(net))], ' not in `net`')
    genes <- genes[which(genes %in% rownames(net))]
    if(length(genes) == 0) {
      return('No gene in `net`')
    }
  }

  multiple <- length(genes) > 1
  if(multiple) {
    y <- Matrix::colMeans(t(scale(t(as.matrix(SingleCellExperiment::logcounts(net)[genes,])))))
  } else {
    v <- SingleCellExperiment::logcounts(net)[genes,]
    y <- (v-mean(v))/(stats::sd(v))
  }
  RR <- do.call(rbind, lapply(seq(1,nRand), function(i) {
    if(multiple) {
      Matrix::colMeans(t(scale(t(as.matrix(SingleCellExperiment::logcounts(net)[sample(rownames(net), length(genes)),])))))
    } else {
      v <- SingleCellExperiment::logcounts(net)[sample(rownames(net), length(genes)),]
      (v-mean(v))/stats::sd(v)
    }
  }))
  RR <- RR[which(apply(RR, 1, function(i) {!any(is.na(i))})),]
  z <- (y-Matrix::colMeans(RR))/MatrixGenerics::colSds(RR)
  p <- 2 * stats::pnorm(-abs(z))

  return(S4Vectors::List(obs=y,
                         z=z,
                         p=p))

}
