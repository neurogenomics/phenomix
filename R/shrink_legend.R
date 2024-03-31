shrink_legend <- function(p,
                          legend.position = "bottom",
                          legend.key.size=0.001,
                          legend.spacing=.0001,
                          legend.text=6,
                          strip.text=legend.text,
                          ncol=5){
    ### make legend smaller
    p + 
        ggplot2::theme(legend.key.size = ggplot2::unit(legend.key.size, "cm"),
                   legend.spacing =ggplot2::unit(legend.spacing,"cm"),
                   legend.text = ggplot2::element_text(size=legend.text),
                   legend.position = legend.position,
                   strip.text = ggplot2::element_text(size=strip.text)
    ) +
        ## make legend only have one column
        ggplot2::guides(fill = ggplot2::guide_legend(ncol = ncol)) 
}