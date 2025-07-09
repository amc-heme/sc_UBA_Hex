#' Compute hexbin information for any supported single-cell object
#'
#' This function computes the hexagonal binning of dimension-reduction coordinates
#' without modifying the original object, enabling on-the-fly computation for any
#' supported object type (e.g., SingleCellExperiment, Seurat).
#'
#' @param object A single-cell object (e.g., SingleCellExperiment, Seurat).
#' @param nbins Integer specifying number of hexagonal bins (default: 80).
#' @param dimension_reduction Character string specifying which reduced-dimension
#'   coordinates to use (default: "UMAP").
#' @param use_dims Integer vector of length two specifying which dimensions
#'   to use from the reduction (default: c(1, 2)).
#'
#' @return A list with elements `cID` (cell assignements) and `hexbin.matrix`.
#' @export
#' @importFrom hexbin hexbin hcell2xy
#' @import SCUBA
get_hexbin_info <- function(object,
                            nbins = 80,
                            dimension_reduction = "UMAP",
                            use_dims = c(1, 2)) {
    # dispatch to method
    UseMethod("get_hexbin_info")
}

#' @export
get_hexbin_info.SingleCellExperiment <- function(object,
                                                 nbins = 80,
                                                 dimension_reduction = "UMAP",
                                                 use_dims = c(1, 2)) {
    if (!dimension_reduction %in% reducedDimNames(object)) {
        stop("Specified dimension reduction not found in SingleCellExperiment object.")
    }
    dr <- reducedDim(object, dimension_reduction)
    .make_hexbin_helper(dr, nbins, use_dims)
}

#' @export
get_hexbin_info.Seurat <- function(object,
                                   nbins = 80,
                                   dimension_reduction = "UMAP",
                                   use_dims = c(1, 2)) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Seurat package required for get_hexbin_info on Seurat objects.")
    }
    # extract embeddings
    if (!dimension_reduction %in% names(object@reductions)) {
        stop(sprintf("Specified dimension reduction '%s' not found in Seurat object.", dimension_reduction))
    }
    dr <- Seurat::Embeddings(object, reduction = dimension_reduction)
    .make_hexbin_helper(dr, nbins, use_dims)
}
