#' Run \pkg{variancePartition} 
#' 
#' Run \link[variancePartition]{fitExtractVarPartModel} 
#' to evaluate the effect of each metadata variable on your data.
#' 
#' @param obj Matrix or \pkg{Seurat} object to run  
#' \link[variancePartition]{fitExtractVarPartModel} on. 
#' @param metadata \link[data.frame]{data.frame} containing metadata
#' to use in \code{form}. If \code{obj} is a \pkg{Seurat} object,
#' this can be left \code{NULL} and metadata 
#' will be extracted automatically.
#' @param is_opengwas Whether the data comes from
#'  \href{https://gwas.mrcieu.ac.uk/}{OpenGWAS}. 
#' If \code{TRUE}, a predefined formula and metadata 
#' processing procedure will be used.
#' @param nCores Number of cores to use in parallel. 
#' @param ... Additional arguments passed to 
#' \link[variancePartition]{fitExtractVarPartModel}.
#' @inheritParams variancePartition::fitExtractVarPartModel
#' 
#' @importFrom BiocParallel register SnowParam
#' @importFrom methods is
#' @importFrom variancePartition fitExtractVarPartModel
#' 
#' @export
run_variancePartition <- function(obj,
                                  metadata=NULL,
                                  nCores=NULL,
                                  form=NULL,
                                  is_opengwas=FALSE,
                                  show_plot=TRUE,
                                  ...){  
    
    if(is.null(form) & (!is_opengwas)){
        stop("Must provide form or set is_opengwas=TRUE (when applicable).")
    }
    #### Register cores ####
    if(is.null(nCores)) nCores <- assign_cores(worker_cores = nCores)$worker_cores
    BiocParallel::register(BiocParallel::SnowParam(nCores)) 
    #### Extract matrix ####
    mat <- extract_matrix(obj = obj)
    #### Extract metadata ####
    if(methods::is(obj,"Seurat")){
        metadata <- extract_metadata(obj = obj) 
    }
    #### Prepare formula #####
    if(is_opengwas){
        form <- ~ log10(N) + log10(nsnp) +   
            (1|population) + (1|build_inferred) +
            (1|priority) + (1|author) + (1|pmid) + 
            (1|category) + (1|subcategory) + 
            (1|year)
        messager("Using predefined formula for OpenGWAS metadata:\n",
                 c(form))
        #### Prepare metadata ##### 
        metadata <-  metadata[colnames(mat),]
        metadata$year <- factor(metadata$year)
        metadata$priority <- factor(metadata$priority, 
                                   levels = sort(unique(metadata$priority)), 
                                   ordered = TRUE) 
        metadata[is.na(metadata)] <- "NA" # Make NA its own category
    }  
    if(is.null(form)){
        stop("Must provide formula to form argument.")
    }  
    #### Run VP ####
    varPart <- variancePartition::fitExtractVarPartModel(exprObj = mat,
                                                         formula = form,
                                                         data = metadata,
                                                         ...)
    
    if(show_plot){
        
    }
    
    
    
    return(varPart)
}
