#' Run \pkg{variancePartition}
#'
#' Run \link[variancePartition]{fitExtractVarPartModel}
#' to evaluate the effect of each metadata variable on your data.
#'
#' @param obj Matrix or \pkg{Seurat} object to run
#' \link[variancePartition]{fitExtractVarPartModel} on.
#' @param metadata \code{data.frame} containing metadata
#' to use in \code{formula}. If \code{obj} is a \pkg{Seurat} object,
#' this can be left \code{NULL} and metadata
#' will be extracted automatically.
#' @param is_opengwas Whether the data comes from
#'  \href{https://gwas.mrcieu.ac.uk/}{OpenGWAS}.
#' If \code{TRUE}, a predefined formula and metadata
#' processing procedure will be used.
#' @param show_plot Print plot. 
#' @param ... Additional arguments passed to
#' \link[variancePartition]{fitExtractVarPartModel}.
#' @inheritParams variancePartition::fitExtractVarPartModel
#' @inheritParams BiocParallel::MulticoreParam
#' 
#' @returns \code{varPart} results.
#' 
#' @export
#' @importFrom BiocParallel register SnowParam
#' @importFrom methods is
#' @importFrom variancePartition fitExtractVarPartModel
run_variancePartition <- function(obj,
                                  metadata = NULL,
                                  workers = NULL,
                                  formula = NULL,
                                  is_opengwas = FALSE,
                                  show_plot = TRUE,
                                  ...) {
    if (is.null(formula) & (!is_opengwas)) {
        stopper("Must provide formula or set is_opengwas=TRUE (when applicable).")
    }
    #### Register cores #### 
    BPPARAM <- assign_cores(workers = workers)
    #### Extract matrix ####
    mat <- scKirby::get_x(obj = obj)
    #### Extract metadata ####
    if (methods::is(obj, "Seurat")) {
        metadata <-  scKirby::get_obs(obj = obj)
    }
    #### Prepare formula #####
    if (is_opengwas) {
        formula <- ~ log10(N) + log10(nsnp) +
            (1 | population) + (1 | build_inferred) +
            (1 | priority) + (1 | author) + (1 | pmid) +
            (1 | category) + (1 | subcategory) +
            (1 | year)
        messager(
            "Using predefined formula for OpenGWAS metadata:\n",
            c(formula)
        )
        #### Prepare metadata #####
        metadata <- metadata[colnames(mat), ]
        metadata$year <- factor(metadata$year)
        metadata$priority <- factor(metadata$priority,
            levels = sort(unique(metadata$priority)),
            ordered = TRUE
        )
        metadata[is.na(metadata)] <- "NA" # Make NA its own category
    }
    if (is.null(formula)) {
        stopper("Must provide formula argument.")
    }
    #### Run VP ####
    varPart <- variancePartition::fitExtractVarPartModel(
        exprObj = mat,
        formula = formula,
        data = metadata,
        BPPARAM = BPPARAM,
        ...
    )
    #### Plot ####
    if (isTRUE(show_plot)) {
        out_merged <- plot_variancePartition(varPart = varPart, 
                                             plot_PercentBars = TRUE, 
                                             plot_VarPart = TRUE)
    } 
    return(varPart)
}
