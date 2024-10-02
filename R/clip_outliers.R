#' Clip outliers
#' 
#' Clip outliers from a ggplot object using Mahalanobis distance.
#' @param gg A ggplot or patchwork object.
#' @param max_dist Maximum Mahalanobis distance to consider as a non-outlier.
#' @inheritParams patchwork::plot_layout
#' @inheritDotParams patchwork::wrap_plots
#' @return A ggplot or patchwork object with outliers clipped.
#' @export
#' @examples
#' gg <- ggplot2::ggplot(mtcars,ggplot2::aes(x=mpg,y=disp))+ggplot2::geom_point()
#' clip_outliers(gg)
clip_outliers <- function(gg,
                          max_dist=NULL,
                          guides = "collect", 
                          axes="collect",
                          axis_titles="collect",
                          ...){  
    
    if(isFALSE(max_dist)) return(gg)
    
    clip_outliers_i <- function(gg, max_dist){
        is.outlier <- mahal.dist <- NULL;
        
        df <- ggplot2::ggplot_build(gg)$data[[1]][,c("x","y")] 
        outlier_md <- rstatix::mahalanobis_distance(data.table::data.table(df))
        if(!is.null(max_dist)){
            outlier_md[,is.outlier:=mahal.dist>max_dist]
        }  
        messager(formatC(sum(outlier_md$is.outlier),big.mark=","),
                 "outliers identified.")
        messager("max_dist for non-outlier:",
                 max(outlier_md[is.outlier==FALSE,]$mahal.dist))
        # ggplot(aes(x,y),data = outlier_md)+ 
        #   geom_point(aes(color=factor(is.outlier)),size=5)    
        gg +
            ggplot2::lims(x=c(min(outlier_md[is.outlier==FALSE]$x),
                              max(outlier_md[is.outlier==FALSE]$x)
            ),
            y=c(min(outlier_md[is.outlier==FALSE]$y),
                max(outlier_md[is.outlier==FALSE]$y)
            )
            )
    }
    if(methods::is(gg,"patchwork")){
        lapply(seq(length(gg)), function(i){
            clip_outliers_i(gg[[i]],max_dist)
        }) |> 
            patchwork::wrap_plots(...) +
            patchwork::plot_layout(guides = guides, 
                                   axes = axes, 
                                   axis_titles =axis_titles)
    } else{
        clip_outliers_i(gg,max_dist)
    }
}