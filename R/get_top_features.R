

get_top_features <- function(seurat,
                             reduction=NULL, 
                             n_features=10,
                             n_features_plot=20, 
                             factors_plot=1:20,
                             invert_vars=F,
                             plot_hist=T,
                             show_plot=T,
                             plot_title=NULL){
    if(is.null(reduction)) reduction <- Seurat::Reductions(seurat)[1]
    message("+ Using reduction: ",reduction)
    loadings <- seurat@reductions[[reduction]]@feature.loadings 
    top_genes <- reshape2::melt(loadings) %>%  
        dplyr::rename(gene=Var1, factor=Var2, loading=value) %>% 
        dplyr::group_by(factor) %>%
        dplyr::slice_max(order_by = abs(loading), n = n_features)  %>%
        data.table::data.table()
    
    if(show_plot){
        gp <- plot_top_phenos(top_genes,  
                              factors = factors_plot,
                              invert_vars=invert_vars, 
                              x = "gene",
                              y = "loading", 
                              fil = "factor", 
                              title = plot_title) 
    }
    return(top_genes)
}