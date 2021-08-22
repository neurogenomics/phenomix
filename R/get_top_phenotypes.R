#' Plot the top phenotypes
#' 
#' Create a radar chart of the phenotypes with the 
#' highest loadings for each reduction factor.
#' 
#' @param obj \pkg{Seurat} object or dimensionality reduction object. 
#' @param metadata Phenotype metadata. 
#' Not needed if \code{obj} is a \pkg{Seurat} object.
#' @param n_phenotypes Number of top phenotypes per reduction factor to select. 
#' @param verbose Print messages.
#' @inheritParams plot_top
#' @inheritParams extract_embeddings
#' 
#' @return \code{data.table} of top phenotypes 
#' @export
#' @importFrom Seurat Reductions
#' @importFrom reshape2 melt
#' @importFrom dplyr %>% rename group_by slice_max
#' @importFrom data.table data.table
#' @examples 
#' data("DEGAS_seurat")
#' top_phenos <- get_top_phenotypes(obj=DEGAS_seurat)
get_top_phenotypes <- function(obj,
                               metadata=NULL,
                               reduction=NULL, 
                               n_phenotypes=3, 
                               factors_plot=1:10,
                               invert_vars=FALSE, 
                               show_plot=TRUE,
                               title=NULL,
                               x="phenotype",
                               y="loading",
                               verbose=TRUE){
    embeddings <- extract_embeddings(obj = obj, 
                                     reduction = reduction,
                                     verbose = verbose)
    if(is.null(metadata)) metadata <- extract_metadata(obj = obj)
    top_phenos <- reshape2::melt(embeddings) %>% 
        merge(metadata, 
              by.x = "Var1", by.y=0, all.x = TRUE) %>% 
        dplyr::rename(phenotype=Var1, factor=Var2, loading=value) %>% 
        dplyr::group_by(factor) %>%
        dplyr::slice_max(order_by = abs(loading), n = n_phenotypes)  %>%
        data.table::data.table()
    
    if(show_plot){
        gp <- plot_top(top_data = top_phenos,
                       x=x,
                       y=y,
                       factors = factors_plot,
                       invert_vars = invert_vars,
                       title = title) 
    }
    return(top_phenos)
}