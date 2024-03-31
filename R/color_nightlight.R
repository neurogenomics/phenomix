color_nightlight <- function(n = 40){
    ggplot2::scale_color_gradientn(
        colors=pals::gnuplot(n = n)[-(seq(as.integer(n*0.2)))])
}