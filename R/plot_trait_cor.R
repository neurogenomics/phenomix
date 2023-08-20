#' Plot top trait-trait correlations
#'
#' Takes the output of \link[phenomix]{find_neighbors} as input to \code{knn}. 
#' @param knn A melted similarity matrix or K-nearest neighbor graph.
#' @param top_n The number of top correlations to plot.
#' @param non_self Remove trait2 that are present in trait1.
#' @param group_var Column in \code{knn} data to group by 
#' when computing \code{top_n} traits.
#' @param type Plot type ("heat" or "bar").
#' @param show_plot Whether to print the plot.
#'
#' @return \code{ggplot} object.
#'
#' @export 
#' @import data.table
#' @importFrom utils head 
#' @examples
#' obj <- get_HPO()[seq(100),]
#' knn <- find_neighbors(obj = obj,
#'                       var1_search = "cardio",
#'                       label_col = "HPO_label")
#' gg_cor <- plot_trait_cor(knn=knn,  top_n = 3)
plot_trait_cor <- function(knn,
                           top_n = 3,
                           non_self = TRUE,
                           group_var = "trait1",
                           show_plot = TRUE,
                           type = "heat") {

    requireNamespace("ggplot2")
    if(type=="heat") requireNamespace("heatmaply")
    trait2 <- trait1 <- NULL; 
    
    if (isTRUE(non_self)) {
        plot_dat <- subset(knn, !trait2 %in% trait1)
    } else {
        plot_dat <- knn
    }
    if (!is.null(top_n)) {
        if(!is.null(group_var)){
            plot_dat <- data.table::setDT(plot_dat)[
                order(abs(get("similarity")),decreasing = TRUE),
                utils::head(.SD, top_n),
                by = get(group_var)]
        } else {
            plot_dat <- plot_dat[seq(1, top_n)]
        } 
    } 
    if(type=="bar"){
        gg_cor <- plot_trait_cor_barplot(plot_dat = plot_dat)    
    } else if(type=="heat"){ 
        mat <- melt_to_mat(dat = knn)
        ### Subset original knn to just the top selected traits ####
        mat <- mat[unique(plot_dat$trait1),
                   unique(plot_dat$trait2), drop=FALSE] 
        #### Draw dendrogram only on axes with >1 sample ####
        dendrogram <- if(sum(dim(mat)>1)==2){
            "both"
        } else if (nrow(mat)==1){"column"} else {
            "row"
        }
        gg_cor <- heatmaply::heatmaply(x = mat, 
                                       dendrogram = dendrogram)
    }
    if (isTRUE(show_plot)) print(gg_cor)
    return(gg_cor)
}

