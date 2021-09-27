
#' Identify marker features for each phenotype
#'
#' Identify marker features for a set of phenotypes
#' by taking the intersection between
#' the top \code{n_features} in two assays.
#'
#' @param seurat \code{Seurat} object.
#' @param terms A list of substrings to search for in column names.
#' @param assay1 First assay name.
#' @param assay1_slot Slot to use from \code{assay1}.
#' @param assay2 Second assay name.
#' @param assay2_slot Slot to use from \code{assay2}.
#' @param n_features Number of features to select per assay.
get_phenotype_markers <- function(seurat,
                                  terms = NULL,
                                  assay1 = "MAGMA",
                                  assay1_slot = "scale.data",
                                  assay2 = "specificity",
                                  assay2_slot = "counts",
                                  n_features = 100) {
    trait_set <- grep(paste(terms, collapse = "|"), colnames(seurat),
        ignore.case = TRUE, value = TRUE
    )
    print(paste(length(trait_set), "traits identified."))
    if (!is.null(assay1)) {
        expressed_genes <- head(sort(rowMeans(Seurat::GetAssayData(seurat, assay = assay1, slot = assay1_slot)[, trait_set]),
            decreasing = TRUE
        ), n_features)
    } else {
        expressed_genes <- NULL
    }
    if (!is.null(assay2)) {
        specific_genes <- head(sort(rowMeans(Seurat::GetAssayData(seurat, assay = assay2, slot = assay2_slot)),
            decreasing = TRUE
        ), n_features)
    } else {
        specific_genes <- NULL
    }


    #### Return marker genes ####
    if (sum(!is.null(assay1), !is.null(assay2)) == 2) {
        marker_genes <- intersect(names(expressed_genes), names(specific_genes))
    } else {
        marker_genes <- if (!is.null(assay1)) expressed_genes else specific_genes
    }
    print(paste(length(marker_genes), "marker genes identified."))
    return(marker_genes)
}
