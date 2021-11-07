#' Plot top trait-trait correlations
#'
#' Takes the output of \link[phenomix]{find_neighbors} as input to \code{knn}.
#'
#' @param knn A melted similarity matrix or K-nearest neighbor graph.
#' @param top_n The number of top correlations to plot.
#' @param non_self Remove trait2 that are present in trait1.
#' @param show_plot Whethe to print the plot.
#'
#' @return \code{ggplot} object.
#'
#' @export
#' @import ggplot2
#' @importFrom data.table setDT
#' @importFrom utils head 
plot_trait_cor <- function(knn,
                           top_n = 10,
                           non_self = TRUE,
                           group_var = "trait1",
                           show_plot = TRUE,
                           type = "heat") {

    trait2 <- trait1 <- similarity <- NULL;
    
    if(type=="heat") requireNamespace("heatmaply")
    if (non_self) {
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
                   unique(plot_dat$trait2)]
        gg_cor <- heatmaply::heatmaply(x = mat)
    }
    if (show_plot) print(gg_cor)
    return(gg_cor)
}

