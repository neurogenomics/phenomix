#' Get the Human Phenotype Ontology
#'
#' Import the Human Phenotype Ontology (HPO) as a Seurat object
#' with precalculated variable genes, specificity stored as an assay,
#' and PCA/UMAP reductions.
#'
#' @source
#' \code{
#' save_path <-  file.path("/Desktop/phenome_decomposition",
#'                         "raw_data/HPO/HPO_seurat.rds")
#' obj <- prepare_HPO(save_path = save_path)
#' piggyback::pb_upload(file = save_path,
#'                      repo = "neurogenomics/phenomix",
#'                      overwrite = TRUE)
#' }
#' @return \pkg{Seurat} object
#' @export
get_HPO <- function() {
    tmp <- get_data(fname = "HPO_seurat.rds")
    obj <- readRDS(tmp)
    return(obj)
}
