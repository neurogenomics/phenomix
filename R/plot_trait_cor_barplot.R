plot_trait_cor_barplot <- function(plot_dat){
    ggplot(
        plot_dat,
        aes(x = trait1, y = similarity, fill = similarity)
    ) +
        geom_bar(stat = "identity") +
        facet_grid(facets = trait2 ~ .) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.background = element_rect(fill = "white"),
            strip.text.y = element_text(angle = 0)
        )
}