#' Plot enrichment of terms in a DAG ontology
#' 
#' Plot enrichment of terms in a DAG ontology.
#' @inheritParams ggrepel::geom_label_repel
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
                            sort_by=c("log2_fold_enrichment"=-1),
                            cluster_col="seurat_clusters",
                            id_col=NULL,
                            label_col="dag_enrich.name_wrap",
                            cluster_id_alpha=.9,
                            cluster_id_size=3,
                            color_col=cluster_col,
                            replace_char=list("."=":",
                                              "_"=":"),
                            p_threshold=.05,
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
                            max.overlaps = 100,
                            na.value = "grey50",
                            force_new=FALSE,
                            return_obj=FALSE,
                            show_plot=TRUE,
                            ...){
    term <- depth <- dag_enrich.name <- dag_enrich.name_wrap <- NULL;
    if("dag_enrich.name" %in% colnames(obj@meta.data) && 
       isFALSE(force_new)){
        messager("Using existing enrichment results in obj.",
                 "To rerun, set force_new=TRUE.")
        #### Reconstruct cluster_enrich results from metadata ####
        dag_cols <- grep("^dag_enrich\\.",colnames(obj@meta.data), value = TRUE)
        cluster_enrich <- data.table::data.table(
            unique(obj@meta.data[,c(cluster_col,dag_cols)])
        )
    } else {
        #### Remove any existing "dag_enrich." columns ####
        dag_cols <- grep("^dag_enrich\\.",colnames(obj@meta.data), value = TRUE)
        if(length(dag_cols)>0) obj@meta.data[,dag_cols] <- NULL
        #### Run enrichment ####
        cluster_enrich <- run_dag_enrich_obj(obj=obj, 
                                             ont=ont,
                                             id_col=id_col,
                                             cluster_col=cluster_col,
                                             min_hits=min_hits,
                                             min_offspring=min_offspring,
                                             replace_char=replace_char,
                                             p_threshold=p_threshold,
                                             q_threshold=q_threshold,
                                             sort_by=sort_by) 
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
        cluster_enrich <- cluster_enrich[,.SD[1],by=c("group")]
        names(cluster_enrich) <- paste0("dag_enrich.",names(cluster_enrich))
        cluster_enrich[,dag_enrich.name_wrap:=stringr::str_wrap(
            dag_enrich.name,
            width = label_width)]
        obj@meta.data$id <- colnames(obj)
        meta <- merge(obj@meta.data, 
                      cluster_enrich,
                      by.x="seurat_clusters",
                      by.y="dag_enrich.group", 
                      all.x=TRUE)
        rownames(meta) <- meta$id
        obj@meta.data <- meta[colnames(obj),]
    }  
    obj <- add_cluster_colors(obj)
    cluster_colors <- add_cluster_colors(obj, return_dict = TRUE)
    dr <- Seurat::DimPlot(obj,
                          alpha = point_alpha,
                          pt.size = point_size,
                          na.value = na.value,
                          label = TRUE,
                          label.color = ggplot2::alpha("white",cluster_id_alpha),
                          label.size = cluster_id_size,
                          group.by = color_col, 
                          reduction = reduction) +
        Seurat::NoLegend() +
        ggplot2::scale_color_manual(values = cluster_colors) +
        theme_nightlight() 
    if(isTRUE(add_labels)){
        dark_theme <- dr[[1]]$theme$plot.background$fill=="black" 
        cluster_enrich2 <- merge(dr[[1]]$layers[[2]]$data, 
                                 cluster_enrich, 
                                 by.x=cluster_col, 
                                 by.y="dag_enrich.group")
        # Seurat::LabelClusters ### <-- Not very flexible and has uninformative errors.
        dr <- dr +
         ggrepel::geom_label_repel(
             data = cluster_enrich2, 
             ggplot2::aes(x=!!ggplot2::sym(names(cluster_enrich2)[[2]]),
                          y=!!ggplot2::sym(names(cluster_enrich2)[[3]]),
                          label=!!ggplot2::sym(label_col),
                          fill=!!ggplot2::sym(color_col)
                          ),
             color=if(dark_theme) "white" else "black",
             arrow = ggplot2::arrow(length = ggplot2::unit(0.01, "npc"),
                                    type="closed"), 
             size = label_size,
             box.padding = box.padding,
             point.padding = ggplot2::unit(0.02, "npc"),
             max.overlaps = max.overlaps,
             min.segment.length = 0,
             ...
             ) +
            ggplot2::scale_fill_manual(
                values = ggplot2::alpha(cluster_colors,label_alpha)
            ) 
    }
    #### Show #### 
    if(isTRUE(show_plot)) methods::show(dr)
    #### Return ####
    out <- list(plot=dr,
                data=cluster_enrich)
    if(isTRUE(return_obj)) out$obj <- obj
    return(out)
}
