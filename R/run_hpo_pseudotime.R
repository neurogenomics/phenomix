run_hpo_pseudotime <- function(obj,
                               disease_ids,
                               dt_genes,
                               learn_graph_control,
                               title=NULL,
                               subtitle=NULL,
                               color_by_symptoms=TRUE,
                               symptom_color="red",
                               bg_colors=c("#5000ff","white"),
                               trajectory_graph_color=ggplot2::alpha("white",.8),
                               point_alpha=.5,
                               ...){
    requireNamespace("monocle3")
    requireNamespace("SeuratWrappers")
    library(dplyr)
    gene_symbol <- disease_id <- NULL;
    
    disease_ids <- unique(disease_ids)
    messager("Running pseudotime analysis for:", 
             paste(disease_ids, collapse = "; "))
    p2d <- dt_genes[disease_id %in% disease_ids,
                    list(count=data.table::uniqueN(gene_symbol)),
                    by=c("hpo_id","disease_id")]#|>
    # data.table::dcast.data.table(formula = hpo_id~disease_id,
    #                              value.var = "count")
    if(nrow(p2d)==0) {
        stopper("No hpo_ids found for the given disease_ids.")
    } else {
        messager("Found", nrow(p2d), "hpo_ids for the given disease_ids.")
    }
    # p2d <- p2d[colnames(ref),]
    # ## Using annotation overlap (258276 p2d pairs)
    # nrow(unique(dt_annot[,c("hpo_id","disease_id")]))
    # p2d <- dt_annot[,count:=1]|>
    #     data.table::dcast.data.table(formula = hpo_id~disease_id,
    #                                  value.var = "count",
    #                                  )
    #### Merge with rest of annotations ####
    ## make Seurat way too slow having this many columns...
    # dt_annot_melt <- data.table::merge.data.table(dt_annot_melt,
    #                                               p2d,
    #                                               by.x = "id",
    #                                               by.y = "hpo_id",
    #                                               all.x = TRUE)  
    hpo_ids <- intersect(unique(p2d$hpo_id),
                         colnames(obj))
    if(length(hpo_ids)==0) {
        stopper("No hpo_ids overlapping with samples (colnames) in obj.")
    } else {
        messager("Found", length(hpo_ids), "matching hpo_ids within obj.")
    }
    cds <- suppressWarnings(
        SeuratWrappers::as.cell_data_set(obj)
    )
    
    k <- min(ncol(cds[,hpo_ids]),
             max(as.integer(cds@colData$seurat_clusters))) -2
    cds_sub <- monocle3::cluster_cells(
        cds = cds[,hpo_ids],
        k = k)
    cds_sub <- monocle3::learn_graph(cds_sub, 
                                     learn_graph_control=learn_graph_control)
    cds_sub <- monocle3::order_cells(cds_sub, root_cells = colnames(cds_sub))
    monocle3::principal_graph(cds) <- monocle3::principal_graph(cds_sub)
    monocle3::principal_graph_aux(cds) <- monocle3::principal_graph_aux(cds_sub) 
    #### Contellation plot (dark) ####
    if(isTRUE(color_by_symptoms)){
        cds@colData$is_symptom <-  as.factor(colnames(cds) %in% hpo_ids)
        plt <- monocle3::plot_cells(cds,
                                    color_cells_by = "is_symptom",
                                    label_cell_groups = FALSE,
                                    trajectory_graph_color = trajectory_graph_color,
                                    alpha = point_alpha,
                                    ...
        ) +
            ggplot2::scale_colour_manual(values = bg_colors) +
            # ggplot2::geom_density_2d_filled(data = . %>% dplyr::filter(is_symptom==TRUE), 
            #                                 contour_var = "ndensity")+
            # ggplot2::scale_fill_manual(
            #     values = c("transparent",
            #                ggplot2::alpha(pals::gnuplot(15),.7))) +
            ggplot2::geom_point(data = . %>% dplyr::filter(is_symptom==TRUE),
                                color=symptom_color,
                                alpha=point_alpha,
                                size=1) +
            theme_nightlight() +
            ggplot2::theme(legend.position = "none")
    } else{
        plt <- monocle3::plot_cells(cds,
                                    ...
        ) 
    }
    plt <- plt + ggplot2::labs(title = title, subtitle = subtitle)
    return(
        list(data=cds,
             plot=plt)
    )
}