#' OpenGWAS metadata
#' 
#' Prepare metadata for GWAS/QTL datasets from the OpenGWAS database.
#' @param meta Metadata for OpenGWAS datasets.
#' @param ids A list of OpenGWAS dataset IDs.
#' @param host phenomix database host URL.
#' @param verbose Print messages.
#' @returns data.table
#' 
#' @export
#' @import data.table
#' @examples
#' meta <- opengwas_meta() 
opengwas_meta <- function(meta = phenomix::OpenGWAS,
                          ids = NULL,
                          host = paste0(
                              "https://phenomix.dsi.ic.ac.uk/",
                              "MAGMA_Files_Public/data/GWAS_munged/"),
                          verbose = TRUE){
    # devoptera::args2vars(opengwas_meta)
    url <- id <- magma_annot <- magma_out <- magma_raw <- NULL;
    
    if(!is.null(ids)) meta <- meta[id %in% ids,]
    if(!"url" %in% names(meta)){
        meta[,url:=paste(host,id,paste0(id,".tsv.bgz"),sep="/")]   
        add_magma_url <- function(host,i,suffix){
            paste0(host,i,"/MAGMA_Files/",
                   paste0(i,".tsv.bgz.35UP.10DOWN/"),
                   paste0(i,".tsv.bgz.35UP.10DOWN.genes",suffix))
        } 
        meta[,magma_annot:=add_magma_url(host,id,suffix=".annot")]   
        meta[,magma_out:=add_magma_url(host,id,suffix=".out")]   
        meta[,magma_raw:=add_magma_url(host,id,suffix=".raw")]    
    }
    messager("Querying",nrow(meta),"dataset ids.",v=verbose)
    # meta <- meta[RCurl::url.exists(url),]
    return(meta)
}