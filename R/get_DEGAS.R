#' DeGAs seurat object
#'
#' Decomposition of Genetic Associations (DeGAs) is a methodology first used in
#' Tanigawa et al. 
#' \href{https://www.nature.com/articles/s41467-019-11953-9}{
#' (2019, Nature Communications)} to reduce 2,138 GWAS to 100 components. 
#' DeGAs is essentially the application of 
#' Truncated Singular Value Decomposition (TSVD) to pre-filtered
#'  GWAS summary statistics across many traits.
#' Here, we have reprocessed this data and stored it as a Seurat object with 
#'  trait embeddings and gene loadings (stored in \code{reductions} slot)
#'   for the gene-level version of the decomposition
#'  ("contributionGene"), where SNPs were aggregated to
#'   gene-level by proximity. Trait metadata is stored in the 
#'   \code{meta.data} slot. We also reconstructed an approximation of the 
#' original gene x trait matrix via matrix multiplication of the trait
#' embeddings matrix (U) and the gene loadings matrix (V) and stored it as the
#' assay "genes".
#' 
#' Specifically, we used the
#'  \emph{dev_allNonMHC_z_center_p0001_100PCs_20180129.npz} file found on the 
#'  \href{https://biobankengine.stanford.edu/degas#download}{
#'  Global Biobank Engine}.
#' 
#' @source \href{https://biobankengine.stanford.edu/degas#download}{Global Biobank Engine}
#' @source \href{https://github.com/rivas-lab/public-resources/tree/master/uk_biobank/DeGAs}{DeGAs GitHub}
#' @source \href{https://www.nature.com/articles/s41467-019-11953-9}{Original publication}
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
#' #### Use matrix multiplication to reconstruct trait x gene matrix ####
#' M <- (u %*% t(v)) %>% Matrix::t() %>% as("sparseMatrix")
#' degas <- Seurat::CreateSeuratObject(counts = M,
#'                                     assay = "genes", 
#'                                     meta.data = DEGAS_metadata)
#' degas[["contributionGene"]]  <- DEGAS_contributionGene
#'
#' save_path <- file.path(tempdir(),"DEGAS.rds")
#' saveRDS(degas, save_path)
#' piggyback::pb_upload(file = save_path,
#'                      repo = "neurogenomics/phenomix",
#'                      overwrite = TRUE)
#' }
#' @return \pkg{Seurat} object
#' @export
#' @examples 
#' degas <- phenomix::get_DEGAS()
get_DEGAS <- function() {
    tmp <- get_data(fname = "DEGAS.rds")
    obj <- readRDS(tmp)
    return(obj)
}
