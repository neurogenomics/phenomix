#' Plot highlights
#' 
#' Plot highlighted subsets of the Seurat object.
#' @inheritParams patchwork::plot_layout
#' @inheritDotParams Seurat::LabelPoints
#' @export
#' @examples
#' obj <- get_HPO()
#' out <- plot_highlights(obj,
#'                        name_col="HPO_label",
#'                        queries=c("Alzheimer","Parkinson"))
plot_highlights <- function(obj,
                            name_col,
                            queries,
                            joint=TRUE,
                            palettes=list("pals::brewer.blues",
                                          "pals::brewer.reds",
                                          "pals::brewer.purples",
                                          "pals::brewer.greens",
                                          "pals::brewer.oranges"),
                            ignore.case = TRUE,
                            ncol=1,
                            axis_titles = "collect",
                            axes = "collect",
                            guides = "collect",
                            wrap_width=20,
                            show_plot=TRUE,
                            ...){
    force(obj)
    force(queries)
    
    nms <- lapply(stats::setNames(queries,
                                  queries), function(q){
        m <- obj@meta.data[grepl(q,obj@meta.data[[name_col]],
                                 ignore.case = ignore.case),]
        stats::setNames(m[[name_col]],rownames(m))
    })
    ids <- lapply(nms,names)
    colors <- lapply(seq(length(nms)), function(i){
        pfun <- eval(parse(text=palettes[[i]]))
        n <- length(nms[[i]])
        stats::setNames(pfun(n = max(3,n))|>tail(n),
                        names(nms[[i]]))
    })|> `names<-`(names(nms))
    cols.highlight <- lapply(colors,tail,1)
    # cols.highlight <- stringr::str_split(palettes,"[.]",simplify = TRUE)[,2][seq(length(queries))]
    
    if(isTRUE(joint)){
        plt <- (
            Seurat::DimPlot(obj,
                            cells.highlight = ids,
                            cols.highlight  = cols.highlight
            )
        )  |>
            Seurat::LabelPoints(points = unname(unlist(ids)), 
                                labels = stringr::str_wrap(
                                    unname(unlist(nms)),
                                    width = wrap_width
                                ),
                                color=unlist(unname(colors)),
                                repel = TRUE, 
                                ...) 
    } else{
        plt <- lapply(seq(length(nms)), function(i){
            (
                Seurat::DimPlot(obj,
                                cells.highlight = ids[[i]],
                                cols.highlight  = colors[[i]]
                ) +
                    Seurat::NoLegend()
            )  |>
                Seurat::LabelPoints(points = unname(ids[[i]]), 
                                    labels = stringr::str_wrap(
                                        unname(nms[[i]]),
                                        width = wrap_width
                                    ),
                                    color=unname(colors[[i]]),
                                    repel = TRUE,
                                    ...) +
                ggplot2::labs(title=queries[[i]])
        })|>
            patchwork::wrap_plots()+
            patchwork::plot_layout(ncol = ncol,
                                   axis_titles = axis_titles,
                                   axes = axes,
                                   guides = guides)
    }
    #### Darken hues ####
    if(show_plot) methods::show(plt)
    return(plt)
}