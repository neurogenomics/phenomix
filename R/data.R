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
#' 
#' \code{
#' #### Embeddings/loadings 
#' DEGAS_contributionGene <- readRDS("/Desktop/phenome_decomposition/raw_data/DEGAS/contributionGene.rds")
#' 
#' #### Metadata
#' DEGAS_metadata <- read.csv("/Desktop/phenome_decomposition/raw_data/DEGAS/metadata_processed.csv", check.names = FALSE)
#' colnames(DEGAS_metadata) <- gsub(" ","_",colnames(DEGAS_metadata))
#' colnames(DEGAS_metadata) <- gsub("[(]|[)]","",colnames(DEGAS_metadata))
#' rownames(DEGAS_metadata) <- DEGAS_metadata$label_phe_code
#' 
#' empty_mat <- Matrix::sparseMatrix(i = nrow(DEGAS_contributionGene@feature.loadings),
#'                                   j = nrow(DEGAS_contributionGene@cell.embeddings),
#'                                   dimnames = list(rownames(DEGAS_contributionGene@feature.loadings), 
#'                                                   rownames(DEGAS_contributionGene@cell.embeddings) 
#'                                   )
#' )
#' DEGAS_seurat <- Seurat::CreateSeuratObject(counts = empty_mat, assay = "genes", meta.data = DEGAS_metadata)
#' DEGAS_seurat[["contributionGene"]]  <- DEGAS_contributionGene
#' usethis::use_data(DEGAS_seurat, overwrite = TRUE)
#' }
"DEGAS_seurat"


 

#' CellTypeDataset: ctd_BlueLake2018_FrontalCortexOnly 
#' 
#' 
#' @source \href{ https://pubmed.ncbi.nlm.nih.gov/29227469/}{Publication} 
#' @source \href{https://www.nature.com/articles/s41467-019-11953-9}{Publication}
#' 
#' \code{
#' ctd_BlueLake2018_FrontalCortexOnly <- readRDS("/Desktop/model_celltype_conservation/processed_data/EWCE/standardized_CTD/BlueLake2018_FrontalCortexOnly.rds")
#' usethis::use_data(ctd_BlueLake2018_FrontalCortexOnly, overwrite = TRUE)
#' } 
"ctd_BlueLake2018_FrontalCortexOnly"




#' MHC genes
#' 
#' Retrieve genes from the MHC 
#' (Major Histocompatibility Complex) region.
#' 
#' The MHC region tends to be ubiquitously associated with all GWAS,
#' so removing it allows us to focus on non-MHC signals of interest.
#' For further information related to the MHC region, see the sources below.
#' 
#' @source \href{https://www.nature.com/articles/s41467-019-11953-9}{DeGAs}
#'  
#' \code{
#' MHC_genes <- get_mhc_genes()
#' usethis::use_data(MHC_genes, overwrite = TRUE)
#' }
"MHC_genes"

