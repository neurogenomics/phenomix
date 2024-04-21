#' Plot preservation histogram
#' 
#' Plot the correlation between equivalent terms in 
#' high-, mid-, and low-dimensional space.
#' @export
plot_preservation_histo <- function(consist_dt,
                                    lvl="hd",
                                    key=c("hd"="High-dimensional",
                                          "md"="Mid-dimensional",
                                          "ld"="Low-dimensional"),
                                    title=paste(key[lvl],"space"), 
                                    subtitle=paste0(
                                        "(",consist_dt[[paste0("dim_",lvl)]][1],
                                        " dimensions)"),
                                    xlim=c(-1,1)
                                    ){
    id_type <- NULL;
    ## histograms of high-dimensional cor 
    ggplot2::ggplot(consist_dt,
                    ggplot2::aes(x=!!ggplot2::sym(paste0("cor_",lvl)),
                                 fill=id_type, color=id_type))+ 
        ggplot2::geom_density(ggplot2::aes(y=ggplot2::after_stat(scaled)),
                              alpha=.5)+ 
        ## baseline mean
        ggplot2::geom_vline(xintercept = mean(consist_dt[[paste0("baseline_cor_",lvl)]]), 
                            linetype="dashed", color="grey40") + 
        ggplot2::geom_text(x = mean(consist_dt[[paste0("baseline_cor_",lvl)]]), y=.8,
                           label="baseline", angle=90, color="grey40",vjust=1.5,
                           check_overlap = TRUE) + 
        ## group means
        ggplot2::geom_vline(data=consist_dt[,.(mean_cor=mean(get(paste0("cor_",lvl)),
                                                             na.rm=TRUE)),
                                            by=id_type],
                            ggplot2::aes(xintercept=mean_cor, color=id_type))+  
        ggplot2::lims(x=xlim )+
        ggplot2::theme_minimal() +
        ggplot2::labs(title=title,
                      subtitle=subtitle,
                      x="Pearson correlation")
}
