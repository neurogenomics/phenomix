#' Theme nightlight
#' 
#' \link{ggplot2} theme with dark background and white text.
#' @param text Text color.
#' @param background Background color.
#' @inheritDotParams ggplot2::theme
#' @returns ggplot2 theme.
#' @export
theme_nightlight <- function(text="white",
                             background="black",
                             color="",
                             ...){
    
    ggplot2::theme(text=ggplot2::element_text(color=text),
                   plot.background = ggplot2::element_rect(fill=background),
                   panel.background = ggplot2::element_rect(fill=background),
                   axis.text = ggplot2::element_text(color=text),
                   axis.line = ggplot2::element_line(color=text),
                   axis.ticks = ggplot2::element_line(color=text),
                   axis.title = ggplot2::element_text(color=text),
                   panel.grid = ggplot2::element_blank(),
                   
                   legend.background = ggplot2::element_rect(fill=background),
                   legend.text = ggplot2::element_text(color=text),
                   strip.background = ggplot2::element_rect(fill=background, 
                                                            color = text),
                   strip.text = ggplot2::element_text(color=text),
                   ...) 
}