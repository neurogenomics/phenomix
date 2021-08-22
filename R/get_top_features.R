#' Plot the top feature loadings
#' 
#' Plot the top feature loadings for from a given reduction. 
#' 
#' @param obj \pkg{Seurat} object or dimensionality reduction object.
#' @param n_features Number of top features per factor to select.
#' @param n_features_plot Number of top features to plot. 
#' @param verbose Print messages. 
#' @inheritParams plot_top
#' @inheritParams extract_loadings
#' 
#' @return top_features \code{data.table}
#' @export
#' @importFrom reshape2 melt
#' @importFrom dplyr %>% rename group_by slice_max
#' @importFrom data.table data.table
#' 
#' @examples
#' data("DEGAS_seurat")
#' top_features <- get_top_features(obj=DEGAS_seurat)
get_top_features <- function(obj,
                             reduction=NULL, 
                             n_features=3,
                             n_features_plot=20, 
                             factors=seq(1,10),
                             invert_vars=FALSE, 
                             show_plot=TRUE,
                             title=NULL,
                             return_melted=FALSE,
                             verbose=TRUE){
    loadings <- extract_loadings(obj = obj, 
                                 reduction = reduction, 
                                 verbose = verbose)
    top_features <- reshape2::melt(loadings) %>%  
        dplyr::rename(feature=Var1, factor=Var2, loading=value) %>% 
        dplyr::group_by(factor) %>%
        dplyr::slice_max(order_by = abs(loading), n = n_features)  %>%
        data.table::data.table()
    
    if(show_plot){
        gp <- plot_top(top_data = top_features,  
                       factors = factors,
                       invert_vars = invert_vars, 
                       x = "feature",
                       y = "loading", 
                       fill = "factor", 
                       title = title) 
    }
    ##### Return format ####
    if(!return_melted){ 
        top_features <- lapply(colnames(loadings), function(x){
            names(sort(abs(loadings[,x]), decreasing = TRUE)[seq(1,n_features)]) 
        }) %>% `names<-`(colnames(loadings))
    } 
    return(top_features)
}

