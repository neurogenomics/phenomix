#' Plot enrichment of terms in a DAG ontology
#' 
#' Plot enrichment of terms in a DAG ontology.
#' @inheritParams Seurat::LabelClusters
#' @inheritDotParams Seurat::LabelClusters
#' @export
#' @examples
#' ont <- KGExplorer::get_ontology("hp")
#' obj <- get_HPO()
#' out <- plot_dag_enrich(obj, ont)
plot_dag_enrich <- function(obj,
                            ont,
                            min_hits = 3, 
                            min_offspring = 100,
                            min_depth=NULL,
                            exact_depth=NULL,
                            exact_ont_ancestors=FALSE,
                            cluster_col="seurat_clusters",
                            label_col="dag_enrich.name_wrap",
                            color_col=cluster_col,
                            replace_char=list("."=":",
                                              "_"=":"),
                            q_threshold=.05,
                            preferred_palettes="kovesi.cyclic_mygbm_30_95_c78",
                            add_labels=TRUE,
                            reduction=NULL,
                            label_alpha = .5,
                            point_alpha=.4,
                            point_size=.7,
                            label_size = 3, 
                            label_width=20,
                            repel = TRUE, 
                            box = TRUE,
                            box.padding = 1,
                            show_plot=TRUE,
                            ...){
    term <- depth <- NULL;
    
    cluster_enrich <- run_dag_enrich_obj(obj=obj, 
                                         ont=ont,
                                         cluster_col=cluster_col,
                                         min_hits=min_hits,
                                         min_offspring=min_offspring,
                                         replace_char=replace_char,
                                         q_threshold=q_threshold) 
    #### Constrain my ancestors in ontology ####
    if(isTRUE(exact_ont_ancestors)){
        cluster_enrich <- cluster_enrich[term %in% ont@elementMetadata$ancestor]    
    }
    #### Filter on min_depth ####
    if(!is.null(min_depth)){
        cluster_enrich <- cluster_enrich[depth<=min_depth]
    }
    #### Filter on exact_depth ####
    if(!is.null(exact_depth)){
        cluster_enrich <- cluster_enrich[depth==exact_depth]
    }
    #### Take only the top hit per cluster ####
    # cluster_enrich|>data.table::setorderv("z_score",-1)
    cluster_enrich_top <- cluster_enrich[,.SD[1],by=c("group")]
    names(cluster_enrich_top) <- paste0("dag_enrich.",names(cluster_enrich_top))
    obj@meta.data$id <- colnames(obj)
    meta <- merge(obj@meta.data, 
                  cluster_enrich_top,
                  by.x="seurat_clusters",
                  by.y="dag_enrich.group", 
                  all.x=TRUE)
    rownames(meta) <- meta$id
    meta$dag_enrich.name_wrap <- stringr::str_wrap(meta$dag_enrich.name,
                                                   width = label_width)
    obj@meta.data <- meta
    na.value <- "grey50"
    dr <- Seurat::DimPlot(obj,
                          cols = KGExplorer::map_colors(
                              obj@meta.data,
                              as = "dict",
                              columns = label_col,
                              preferred_palettes = preferred_palettes)[[1]],
                          alpha=point_alpha,
                          pt.size = point_size,
                          na.value = na.value,
                          group.by = label_col, 
                          reduction = reduction) +
        Seurat::NoLegend() +
        theme_nightlight()
    if(isTRUE(add_labels)){
        dark_theme <- dr[[1]]$theme$plot.background$fill=="black"
        dr <- Seurat::LabelClusters(
            dr,
            id = label_col,
            box = box,
            fill = ggplot2::alpha(
                setdiff(unique(ggplot2::ggplot_build(dr)$data[[1]]$colour),
                        na.value),
                label_alpha),
            color=if(dark_theme) "white" else "black",
            arrow = ggplot2::arrow(length = ggplot2::unit(0.01, "npc"),
                                   type="closed"),
            repel = repel, 
            size = label_size, 
            box.padding = box.padding,
            ...) 
    }
    if(isTRUE(show_plot)) methods::show(dr)
    return(
        list(obj=obj,
             plot=dr,
             data=cluster_enrich)
    )
}
