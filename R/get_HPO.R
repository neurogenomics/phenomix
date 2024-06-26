#' Get the Human Phenotype Ontology
#'
#' Import the Human Phenotype Ontology (HPO) as a Seurat object
#' with precalculated variable genes, specificity stored as an assay,
#' and PCA/UMAP reductions.
#'
#' @source
#' \code{
#' save_path <- file.path(tempdir(),"prepare_hpo.rds")
#' obj <- prepare_hpo(save_path = save_path)
#' piggyback::pb_upload(file = save_path,
#'                      repo = "neurogenomics/phenomix",
#'                      overwrite = TRUE)
#' }
#' @return \pkg{Seurat} object
#' @export
#' @examples
#' obj <- get_HPO() 
get_HPO <- function() {
    tmp <- get_data(file = "prepare_hpo.rds")
    readRDS(tmp)
}
