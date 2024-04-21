#' Plot ontological velocity
#' 
#' Plot ontological velocity.
#' @inheritParams metR::scale_mag
#' @inheritParams metR::geom_streamline
#' @inheritDotParams Seurat::DimPlot
#' @import data.table
#' @export
#' @examples
#' obj <- get_HPO()
#' ont <- KGExplorer::get_ontology("hp")
#' out <- plot_ontological_velocity(obj,ont)
plot_ontological_velocity <- function(obj,
                                      graph_name = tail(names(obj@graphs),1),
                                      reduction= tail(Seurat::Reductions(obj),1),
                                      ont = KGExplorer::get_ontology("hp"),
                                      id_col=NULL,
                                      iterations=1000, 
                                      max_dist=.2,
                                      k=20,# Matches default from Seurat::FindNeighbors
                                      agg_fun=function(x){mean(x,na.rm=TRUE)},
                                      max_size=1,
                                      linewidth_range=c(0,1),
                                      L = 5,
                                      res = 5, 
                                      show_plot=TRUE,
                                      ...){
    # https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html
    # https://ouyanglab.com/singlecell/clust.html
    # https://www.bioconductor.org/packages/devel/bioc/vignettes/Nebulosa/inst/doc/nebulosa_seurat.html
    # https://eliocamp.github.io/metR/reference/geom_arrow.html
    # brew install boost
    # obj2 <- SeuratWrappers::RunVelocity(obj)
    requireNamespace("ggplot2")
    requireNamespace("metR")
    
    # g <- KGExplorer::ontology_to(ont,"tbl_graph") 
    ## Step 1. Sample K-nearest neighbours for 100 samples 
    ## of the low-dimensional graph embedding (UMAP).
    if(!is.null(id_col)){
        if(id_col=="hp_id"){
            obj$hp_id <- data.table::fcoalesce(
                obj$hp_id,
                ifelse(grepl("^HP",colnames(obj)),
                       colnames(obj),NA))   
        }
        obj <- obj[,!is.na(obj@meta.data[[id_col]])]
        colnames(obj) <- make.unique(map_id_sep(obj$hp_id))
    }   
    slopes <- run_ontological_velocity_grid(obj = obj,
                                            ont = ont, 
                                            graph_name = graph_name,
                                            reduction = reduction,
                                            max_dist = max_dist,
                                            iterations = iterations,
                                            agg_fun = agg_fun)   
    ## Avoids internal bug with metR
    unit <- utils::getFromNamespace(x = "unit", ns = "grid")
    p_arrow <- Seurat::DimPlot(obj, 
                               ...) +
        ggnewscale::new_scale_color() +
        metR::geom_arrow(data=slopes,
                         pivot = .5,# 0 beginning, 1=end
                         # min.mag = 0,
                         # arrow.length=1,
                         alpha=.7,
                         arrow.type = "closed",
                         ggplot2::aes(x=x.id1, y=y.id1,
                                      dx = dx, dy = dy,
                                      # alpha=ont_velocity,
                                      color=ont_velocity
                                      )) +
        metR::scale_mag(max_size = max_size) +
        ggplot2::scale_color_viridis_c()+
        ggplot2::theme_bw()
    
 
    slopes[is.na(dx),dx:=0]
    slopes[is.na(dy),dy:=0]
    data.table::setorderv(slopes,"I")
    p_streamline <- ggplot2::ggplot(slopes,
                    ggplot2::aes(x=x, y=y)
                    ) +
        metR::geom_contour_fill(ggplot2::aes(
            fill = ggplot2::after_stat(level),
            z=ifelse(is.na(ont_velocity),0,ont_velocity)
            )
            ) +
            ggplot2::guides(fill =ggplot2::guide_colorsteps()) +
        metR::geom_streamline(ggplot2::aes(dx = dx, 
                                           dy = dy,
                                           linewidth = ggplot2::after_stat(step),
                                           alpha = ggplot2::after_stat(step),
                                           color = ggplot2::after_stat(sqrt(dx^2 + dy^2))
                                           ),
                              color="white",
                              arrow = NULL,
                              L = L,
                              res = res, 
                              lineend = "round"
                              ) +
            ggplot2::scale_color_viridis_c(option = "plasma", guide = "none") +
            ggplot2::scale_linewidth(range = linewidth_range, guide = "none") +
            ggplot2::theme_bw()
            # metR::scale_fill_divergent_discretised(low = "transparent",mid = "blue",high="red")  
    
    if(isTRUE(show_plot)) methods::show(p_arrow)
    return(list(
        plot=list(
            arrow=p_arrow,
            streamline=p_streamline
        ),
        data=slopes
    ))
}


# Seurat::FeaturePlot(obj, pt.size = 1,
#                     # cols=c("grey","blue","red"),
#                     features  = "ontLvl") +
#     ggplot2::scale_color_viridis_c(option="plasma")

# ggplot2::ggplot() +
# ggplot2::geom_point(data=obj@reductions$umap@cell.embeddings,
#                     ggplot2::aes(
#     x = umap_1,
#     y = umap_2),
#     alpha=.1) +
# ggplot2::geom_curve(data=slopes[abs(dx)>0.001 & abs(dy)>0.001,],
#                       ggplot2::aes(
#     x = umap_1.id1,
#     y = umap_2.id1,
#     xend = umap_1.id2,
#     yend = umap_2.id2,
#     color = abs(shortest_distance),
#     alpha=abs(shortest_distance),
#     linewidth = abs(shortest_distance)
#     ),
#     arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "inches"))
#     )


# lvls <- seq(max(ont@elementMetadata$ontLvl))
# for (l in lvls){
#     ont <- KGExplorer::add_ancestors(ont, 
#                                      prefix = paste0("ancestor",l),
#                                      lvl=l, 
#                                      force_new=TRUE)
# }
# 
# meta <- merge(
#       obj@meta.data,
#       ont@elementMetadata[,],
#       by.x=0,
#       by.y="id",
#       all.x = TRUE
# )
# rownames(meta) <- meta$Row.names
# obj@meta.data <- meta 
# clustree::clustree(obj,
#                    prefix = "ancestor") 
