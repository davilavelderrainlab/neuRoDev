#' Reference signatures correlation
#'
#' @param S The new signatures (a list or a matrix)
#' @param refS Reference signatures (a list or a matrix) with the same rownames as S
#' @param warn A boolean (defaults to TRUE) that specifies whether a warning
#' should be given in case the signatures are too similar within them and thus
#' the results should be interpreted carefully.
#'
#' @return A dataframe that contains the correlation values, the highest correlation
#' value for each cluster/row of S, the best label (based on highest correlation
#' value). In case the refS is a list, the highest value and best label are also
#' specified for each element of the list.
#' @export
#'
#' @examples
#' set.seed(123)
#' S <- FC_signatures(matrix(runif(200,0,10), ncol = 10))
#' rownames(S) <- paste0('Gene-', seq(1, dim(S)[1]))
#' refS <- FC_signatures(matrix(runif(200,0.1,7), ncol = 10))
#' colnames(refS) <- paste0('Reference-', seq(1, dim(refS)[2]))
#' rownames(refS) <- paste0('Gene-', seq(1, dim(refS)[1]))
#' reference_signatures_correlation(S, refS)
reference_signatures_correlation <- function(S,
                                             refS,
                                             warn=TRUE) {

  if (!is.list(S)) S <- list(S)
  if (!is.list(refS)) refS <- list(refS)

  if(warn) {
    mean_sig <- mean(unlist(S) > 1)
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
  }

  if (is.null(names(refS)) & length(refS) > 1) names(refS) <- as.character(seq(1, length(refS)))

  calc_correlation <- function(S, refS) {
    common_genes <- intersect(rownames(S), rownames(refS))
    if (length(common_genes) == 0) return(NULL)

    y1 <- refS[common_genes,]
    i1 <- S[common_genes,]
    stats::cor(y1, i1)
  }

  compute_assignments <- function(S, refS) {
    lapply(seq_along(S), function(i) {
        df_cells <- matrix(nrow = ncol(S[[i]]), ncol = 0)
        rownames(df_cells) <- colnames(S[[i]])
        column_list <- vector("list", length = ncol(S[[i]]))
        names(column_list) <- colnames(S[[i]])

        cor_values <- lapply(seq_along(refS), function(y) {
          c <- calc_correlation(S[[i]], refS[[y]])
          if (is.null(c)) return(NULL)

          t_c <- t(c)
          if (!is.null(names(refS))) {
            colnames(t_c) <- paste(names(refS)[y], colnames(t_c), sep = '.')
          }
          df_cells <- assign('df_cells', cbind(df_cells, t_c), envir = '.GlobalEnv')

          lapply(seq_len(ncol(c)), function(z) {
            res <- c[, z]
            if (!is.null(names(refS))) {
              column_list[[colnames(c)[z]]][[names(refS)[y]]] <<- res
            } else {
              column_list[[colnames(c)[z]]][[y]] <<- res
            }
          })
          return(c)
        })
        S4Vectors::SimpleList('List' = column_list, 'Df' = df_cells)
      })
  }

  all_correlations <- compute_assignments(S, refS)
  dfs <- lapply(all_correlations, function(i) i$Df)
  all_correlations <- lapply(all_correlations, function(i) i$List)

  all_maxima <- lapply(all_correlations, function(i) {
    lapply(i, function(x) {
      lapply(x, function(z) {
        if(all(is.na(z))) {
          v <- 0
          names(v) <- 'NA'
          return(v)
        }
        m <- max(z,
                 na.rm = TRUE)
        names(m) <- names(z)[which.max(z)]
        return(m)
      })
    })
  })

  abs_maxima_label <- lapply(all_maxima, function(i) {
    lapply(i, function(x) {
      res <- names(unlist(x))[which.max(unlist(x))]
      if (length(refS) > 1) {
        num_dots <- length(strsplit(names(x)[which.max(unlist(x))], '.', fixed = TRUE)[[1]])
        res <- paste(strsplit(res, ".", fixed = TRUE)[[1]][seq(num_dots + 1, length(strsplit(res, ".", fixed = TRUE)[[1]]))], collapse = '.')
      }
      return(res)
    })
  })

  all_cts <- names(unlist(all_correlations[[1]][[1]]))

  assignment_df <- lapply(seq_along(abs_maxima_label), function(i) {
    df <- S4Vectors::DataFrame(Cluster = names(abs_maxima_label[[i]]))
    rownames(df) <- df$Cluster
    df["Best.Assignment"] <- unlist(abs_maxima_label[[i]])
    df["Best.Value"] <- unlist(lapply(all_maxima[[i]], function(x) max(unlist(x),
                                                                       na.rm=TRUE)))

    for (n in names(refS)) {
      df[paste0(n, ".Assignment")] <- unlist(lapply(all_maxima[[i]], function(x) names(x[[n]])))
      df[paste0(n, ".Value")] <- unlist(lapply(all_maxima[[i]], function(x) x[[n]]))
    }

    no_na_dfs <- dfs[[i]]
    no_na_dfs[which(is.na(no_na_dfs))] <- 0

    cbind(df, no_na_dfs)
  })

  names(assignment_df) <- names(abs_maxima_label)
  assignment_df <- lapply(assignment_df, function(i) {
    colnames(i) <- gsub(' ', '_', colnames(i))
    colnames(i) <- gsub('-', '.', colnames(i))
    colnames(i) <- gsub('/', '.', colnames(i))
    colnames(i) <- gsub('?', '.', colnames(i), fixed = TRUE)
    colnames(i) <- gsub(',', '.', colnames(i), fixed = TRUE)
    return(i)
  })

  if (length(assignment_df) == 1) {
    assignment_df <- assignment_df[[1]]
  } else {
    assignment_df <- do.call(rbind, assignment_df)
  }

  return(assignment_df)

}
