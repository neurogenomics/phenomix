#' GWAS Atlas \pkg{Seurat} object
#'
#' Contains MAGMA summary statistics ("expression" matrix), GWAS metadata,
#'  dimensionality reductions, clusters, and cluster markers.
#'
#'
#' @family CORCES_2020
#' @source \href{https://atlas.ctglab.nl/}{GWAS Atlas}
#' @examples
#' \dontrun{
#' # path <- file.path("/Desktop/phenome_decomposition/raw_data",
#' #                  "GWAS_Atlas/GWAS_Atlas_all.seurat.rds")
#'
#'
#' #### Metadata ####
#' root <- "/Desktop/phenome_decomposition/raw_data/GWAS_Atlas"
#' meta <- data.table::fread(file.path(root, "gwasATLAS_v20191115_metadata.txt"))
#' meta$Trait_id <- paste(meta$Trait, meta$id, sep = "_")
#' rownames(meta) <- meta$Trait_id
#'
#' magma_genes <- data.table::fread(
#'     file.path(root, "gwasATLAS_v20191115_magma_P.txt.gz"),
#'     header = TRUE
#' )
#' rownames(magma_genes) <- magma_genes$GENE
#' magma_genes$GENE <- NULL
#' mat <- orthogene::aggregate_mapped_genes(
#'     gene_df = magma_genes,
#'     non121_strategy = "mean"
#' )
#' colnames(mat) <- meta$Trait_id
#' #### Invert p-values ####
#' mat <- -log10(mat) # 1 - mat
#' mat[is.na(mat)] <- 0
#' obj <- scNLP::seurat_pipeline(
#'     counts = mat,
#'     meta.data = meta,
#'     nfeatures = 10000,
#'     vars.to.regress = c("N", "Year"), # "Nsnps","SNPh2_z"
#'     add_specificity = TRUE,
#'     assay_name = "MAGMA"
#' )
#' path <- file.path(tempdir(), "GWAS_Atlas_seurat.rds")
#' saveRDS(obj, path)
#' #### piggyback ####
#' piggyback::pb_upload(
#'     file = path,
#'     repo = "neurogenomics/phenomix",
#'     overwrite = TRUE
#' )
#' }
#' @export
get_GWAS_Atlas <- function() {
    tmp <- get_data(file = "GWAS_Atlas_all.seurat.rds")
    dat <- readRDS(tmp)
    return(dat)
}
