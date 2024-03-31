#' Plot feature density
#' 
#' Plot feature density.
#' @param joint_only Only return the joint density plot.
#' @inheritParams patchwork::plot_layout
#' @inheritParams Nebulosa::plot_density
#' @inheritDotParams Nebulosa::plot_density
#' @export
#' @examples
#' obj <- get_HPO()
#' out <- plot_feature_density(obj, features =  c("C5","C7","C8A","C8B"))
plot_feature_density <- function(obj,
                                 features=sample(rownames(obj),3),
                                 alpha=.5,
                                 joint=TRUE,
                                 joint_only=FALSE,
                                 joint_title=if(length(unique(features))>10){
                                     paste0("Density of ",
                                            length(unique(features)),
                                            " features")
                                 },
                                 axis_titles = "collect",
                                 axes = "collect",
                                 guides = NULL,
                                 ...){
    requireNamespace("Nebulosa")
    
    features2 <- intersect(rownames(obj),features)
    if(length(features2)==0){
        stopper("No intersecting features between object and input features.")
    }
    messager(length(features2),"/",length(features),"features found in object.")
    
    plts <- Nebulosa::plot_density(obj,
                                   alpha=alpha,
                                   features =  features2,
                                   joint = joint,
                                   ...) 
    ### Add theme ####
    plts <- lapply(plts,function(p){
        p + theme_nightlight() + 
            color_nightlight()
    })
    if(!is.null(joint_title)){
        plts[[length(plts)]] <- plts[[length(plts)]] + 
            ggplot2::labs(title=joint_title)    
    } 
    subplots <- lapply(seq(length(plts)-1),function(i) plts[[i]])

    if(isTRUE(joint_only)){
        return(plts[[length(plts)]])
    } else{
        patchwork::wrap_plots(subplots)/plts[[length(plts)]] +
            patchwork::plot_layout(axis_titles = axis_titles,
                                   axes = axes,
                                   guides = guides) 
    }
}