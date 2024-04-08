#' Scale color nightlight
#' 
#' Add color gradient with nightlight colors.
#' @export
#' @param reverse Reverse the color gradient.
#' @param alpha Transparency of the colors.
#' @inheritParams pals::gnuplot
#' @inheritParams ggplot2::scale_color_gradientn
#' @inheritDotParams ggplot2::scale_color_gradientn
scale_color_nightlight <- function(n = 3,
                                   colors=pals::gnuplot(n = n+1)[
                                       -(seq(as.integer((n+1)*0.2)))
                                       ],
                                   reverse=FALSE,
                                   alpha=NULL,
                                   ...){
    if(!is.null(alpha)) {
        if(length(alpha)==1){
            colors <- ggplot2::alpha(colors,alpha)    
        } else if(length(alpha)==length(colors)){
            for(i in seq_along(colors)){
                colors[[i]] <- ggplot2::alpha(colors[[i]],alpha[i])
            }
        } else {
            stopper("Length of alpha must be 1 or equal to length of colors.")
        }
    }
    if(isTRUE(reverse)) colors <- rev(colors)
    ggplot2::scale_color_gradientn(colors=colors,
                                   ...)
}