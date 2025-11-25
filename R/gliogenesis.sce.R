#' SingleCellExperiment object
#'
#' Function to automatically download the gliogenesis reference network along
#' with the expression profiles and metadata. If the object was previously
#' downloaded with this function, the function will also directly load the
#' object.
#'
#' @param directory The directory in which to download the object or from which
#' to load it. If NULL (default), it downloads it in the cache of the neuRoDev
#' package.
#'
#' @format A `SingleCellExperiment` object with:
#' \describe{
#'   \item{assays}{A list containing one matrix named "signatures" (genes x clusters),
#'   one matrix named "logcounts" (genes x clusters)}
#'   \item{colData}{A `DataFrame` with clusters-level metadata (stage, subclass,
#'   study etc.)}
#'   \item{metadata}{A list containing one `SingleCellExperiment` object with
#'   subclass-level pseudobulks (`subclass_psb`), one `SingleCellExperiment` object
#'   with stage-level pseudobulks (`stage_psb`), one correlation matrix with
#'   cell type reference signatures (`correlations`), one `List` object with the
#'   network information (`network`)}
#' }
#'
#' @source Obtained by processing the gliogenesis lineage of the single-cell
#' data compendium.
#' @export
gliogenesis.sce <- function(directory=NULL) {
  figshare_url <- "https://figshare.com/ndownloader/files/XXXXX"  # replace

  if(is.null(directory)) {
    directory <- tools::R_user_dir("neuRoDev", "cache")
    if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
  }

  destfile <- file.path(directory, "gliogenesis.sce.rda")

  if (!file.exists(destfile)) {
    message("Downloading the gliogenesis.sce object (~87MB) to:\n", destfile)
    utils::download.file(figshare_url, destfile, mode = "wb")
  }

  obj_name <- load(destfile)
  get(obj_name)
}
