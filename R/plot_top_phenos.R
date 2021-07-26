
plot_top_phenos <- function(top_phenos,  
                            factors=1:20,
                            x="label_phe", 
                            y="loading", 
                            fill="factor",
                            title=NULL,
                            invert_vars=F,
                            show_plot=T){ 
    if(is(factors,"character")) factors <- as.numeric(stringr::str_split(factors,"_", simplify = T)[,2])
    if(invert_vars) {x1=x; fill1=fill; x=fill; fill=x1;} 
    reduction <- stringr::str_split(top_phenos$factor[1],"_")[[1]][[1]] 
    #### Filtering    #### 
    top_phenos <- top_phenos[top_phenos[["factor"]] %in% paste(reduction,factors,sep="_"), ] 
    #### Plot     #### 
    gp <- ggplot(top_phenos, mapping = aes_string(x=x, y=y, fill=fill)) +  
        geom_bar(stat = "identity") +
        coord_polar() +
        theme_minimal() +
        labs(title = title)
    if(show_plot) print(gp)
    return(gp)
}
