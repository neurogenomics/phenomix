

get_top_phenos <- function(seurat,
                           reduction=NULL, 
                           n_phenotypes=10, 
                           factors_plot=1:20,
                           invert_vars=F,
                           plot_hist=T,
                           show_plot=T,
                           plot_title=NULL){
    if(is.null(reduction)) reduction <- Seurat::Reductions(seurat)[1]
    message("+ Using reduction: ",reduction)
    loadings <- seurat@reductions[[reduction]]@cell.embeddings 
    top_phenos <- reshape2::melt(loadings) %>% 
        merge(seurat@meta.data, 
              by.x = "Var1", by.y =0, all.x = T) %>% 
        dplyr::rename(phenotype=Var1, factor=Var2, loading=value) %>% 
        dplyr::group_by(factor) %>%
        dplyr::slice_max(order_by = abs(loading), n = n_phenotypes)  %>%
        data.table::data.table()
    
    if(show_plot){
        gp <- plot_top_phenos(top_phenos = top_phenos,  
                              factors = factors_plot,
                              invert_vars = invert_vars,
                              title = plot_title) 
    }
    return(top_phenos)
}