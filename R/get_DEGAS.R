#' DeGAs seurat object
#'
#' Contains embeddings/loadings for contributionGene,
#' as well as phenotype metadata.
#'
#' \bold{NOTE:} The assay data is intentionally filled with
#' an empty sparse matrix due to size limits.
#'
#' @source \href{https://github.com/rivas-lab/public-resources/tree/master/uk_biobank/DeGAs}{DeGAs GitHub}
#' @source \href{https://www.nature.com/articles/s41467-019-11953-9}{Publication}
#' @source
#' \code{
#' #### Embeddings/loadings
#' DEGAS_contributionGene <- readRDS(
#'     "/Desktop/phenome_decomposition/raw_data/DEGAS/contributionGene.rds")
#'
#' #### Metadata
#' DEGAS_metadata <- read.csv(
#'     "/Desktop/phenome_decomposition/raw_data/DEGAS/metadata_processed.csv",
#'     check.names = FALSE)
#' colnames(DEGAS_metadata) <- gsub(" ","_",colnames(DEGAS_metadata))
#' colnames(DEGAS_metadata) <- gsub("[(]|[)]","",colnames(DEGAS_metadata))
#' rownames(DEGAS_metadata) <- DEGAS_metadata$label_phe_code
#'
#' empty_mat <- Matrix::sparseMatrix(
#'     i = nrow(DEGAS_contributionGene@feature.loadings),
#'     j = nrow(DEGAS_contributionGene@cell.embeddings),
#'     dimnames = list(rownames(DEGAS_contributionGene@feature.loadings),
#'                   rownames(DEGAS_contributionGene@cell.embeddings)
#'                                   )
#' )
#' DEGAS_seurat <- Seurat::CreateSeuratObject(counts = empty_mat,
#'                                            assay = "genes", 
#'                                            meta.data = DEGAS_metadata)
#' DEGAS_seurat[["contributionGene"]]  <- DEGAS_contributionGene
#'
#' save_path <- file.path(tempdir(),"DEGAS.rds")
#' saveRDS(DEGAS_seurat, save_path)
#' piggyback::pb_upload(file = save_path,
#'                      repo = "neurogenomics/phenomix",
#'                      overwrite = TRUE)
#' }
#' @return \pkg{Seurat} object
#' @export
get_DEGAS <- function() {
    tmp <- get_data(fname = "DEGAS.rds")
    obj <- readRDS(tmp)
    return(obj)
}
