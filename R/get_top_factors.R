
get_top_factors <- function(seurat,
                            reduction=NULL,
                            term,
                            search_col="label_phe",
                            n_quantiles=10,
                            select_quantiles=max(n_quantiles),
                            plot_hist=T){
    if(is.null(reduction)) reduction <- Seurat::Reductions(seurat)[1]
    message("+ Using reduction: ",reduction)
    select_cols <- rownames(seurat@meta.data[grepl(term,unname(seurat@meta.data[[search_col]]),
                                                   ignore.case = T),])  
    if(length(select_cols)>0) message("+ ",length(select_cols),"matching phenotypes identified.") else stop("0 matching phenotypes identified.")
    loadings <- seurat@reductions[[reduction]]@cell.embeddings[select_cols,]
    #### Handle multiple phenotypes per term
    if(!is(loadings,"numeric")) loadings <- colMeans(loadings, na.rm = T)
    if(plot_hist) print(hist(loadings, 50, main=paste("Top",reduction,"factors:",term)))
    
    quantiles <- cut(abs(loadings), breaks = n_quantiles, labels = 1:n_quantiles)
    top_factors <- loadings[quantiles %in% select_quantiles]  
    return(top_factors)
}

