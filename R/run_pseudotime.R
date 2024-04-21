run_pseudotime <- function(obj,
                           disease_ids,
                           root_cells,
                           map_root_cells=TRUE,
                           p2d,
                           id_col,
                           learn_graph_control,
                           use_partition=TRUE,
                           color_cells_by = "is_symptom",
                           title=NULL,
                           subtitle=NULL,
                           color_by_symptoms=TRUE,
                           symptom_color="red",
                           bg_colors=c("#5000ff","white"),
                           trajectory_graph_color=ggplot2::alpha("white",.8),
                           point_alpha=.5,
                           add_density=TRUE,
                           density_filled=FALSE,
                           density_adjust=1,
                           density_alpha=point_alpha,
                           ...){
    requireNamespace("monocle3")
    requireNamespace("SeuratWrappers")
    requireNamespace("dplyr")
    disease_id <- mondo_id <- NULL; 
    
    cat("\nRunning pseudotime analysis for:", 
        paste(disease_ids, collapse = "; "),"\n")
    p2d <- p2d[disease_id %in% disease_ids|  
               mondo_id %in% disease_ids,]
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
    hpo_ids <- if(!is.null(id_col) && 
       id_col %in% colnames(obj@meta.data)) {
        obj@meta.data[[id_col]]
    } else {
        colnames(obj)
    }
    hpo_ids <- intersect(unique(p2d$hpo_id),
                         map_id_sep(unique(hpo_ids))
                         )
    if(length(hpo_ids)==0) {
        stopper("No hpo_ids overlapping with samples (colnames) in obj.")
    } else {
        messager("Found", length(hpo_ids), "matching hpo_ids within obj.")
    }
    #### Convert to CDS format ####
    cds <- suppressWarnings(
        SeuratWrappers::as.cell_data_set(obj)
    ) 
    #### Choose root cells ####
    if(length(root_cells)==0) {
        messager("Setting disease_ids as root cells.")
        root_cells <- disease_ids
    }
    root_cells <- as.character(na.omit(root_cells))
    if(length(root_cells)>0){ 
        root_cells <- c(
            colnames(obj)[map_id_sep(colnames(obj)) %in% map_id_sep(root_cells)],
            if(isTRUE(map_root_cells) && 
               !is.null(id_col) &&
               id_col %in% colnames(obj@meta.data)){
                colnames(obj)[map_id_sep(obj@meta.data[[id_col]]) %in% map_id_sep(root_cells)]    
            },
            if(isTRUE(map_root_cells) && 
               "mondo_id" %in% colnames(obj@meta.data)){
                colnames(obj)[map_id_sep(obj@meta.data[["mondo_id"]]) %in% map_id_sep(root_cells)]
            }
        )|> unique() 
    } 
    #### Choose k ####
    k <- min(ncol(cds[,unique(c(hpo_ids,root_cells))]),
             max(as.integer(cds@colData$seurat_clusters), na.rm = TRUE)
             ) - 2
    k <- max(k,2,na.rm = TRUE)
    messager(paste0("Running monocle3::cluster_cells with k=",k))
    cds_sub <- monocle3::cluster_cells(
        cds = cds[,unique(c(hpo_ids,root_cells))],
        k = k)
    messager("Running monocle3::learn_graph.")
    cds_sub <- monocle3::learn_graph(cds_sub, 
                                     use_partition=use_partition,
                                     learn_graph_control=learn_graph_control)
    if(length(root_cells)==0){
        messager("0 root_cells found. Using all cells as roots.")
        root_cells <- colnames(cds_sub)
    } else{
        messager("Using",length(root_cells),"root cells.")
    }
    cds@colData$is_root <- as.factor(colnames(cds) %in% root_cells)
    messager("Running monocle3::order_cells with",ncol(cds_sub),"cells and",
             length(root_cells),"root cells. ") 
    cds_sub <- monocle3::order_cells(cds_sub, 
                                     root_cells = root_cells)
    monocle3::principal_graph(cds) <- monocle3::principal_graph(cds_sub)
    monocle3::principal_graph_aux(cds) <- monocle3::principal_graph_aux(cds_sub) 
    #### Constellation plot (dark) ####
    if(isTRUE(color_by_symptoms)){
        cds@colData$is_symptom <-  as.factor(colnames(cds) %in% hpo_ids)
        plt <- monocle3::plot_cells(cds,
                                    color_cells_by = color_cells_by,
                                    label_cell_groups = FALSE,
                                    trajectory_graph_color = trajectory_graph_color,
                                    alpha = point_alpha,
                                    ...
        ) +
            ggplot2::scale_colour_manual(values = bg_colors) +
            ggplot2::geom_point(data = . %>% dplyr::filter(is_symptom==TRUE),
                                color=symptom_color,
                                alpha=point_alpha,
                                size=1) +
            theme_nightlight() +
            ggplot2::theme(legend.position = "none")
    } else{
        plt <- monocle3::plot_cells(cds,
                                    color_cells_by=color_cells_by,
                                    ...
        ) 
    }
    # plt$layers <- plt$layers[ 1:10]
    # plt
    graph_layers <- seq(3,9)
    if(add_density){ 
        if(density_filled){
            plt <- plt + 
                ggplot2::geom_density2d_filled(data = . %>% dplyr::filter(is_symptom==TRUE),
                             adjust=density_adjust) +
                ggplot2::scale_fill_manual(
                    values = c("transparent",
                               ggplot2::alpha(pals::gnuplot(15),density_alpha))
                ) 
        } else {
            plt <- plt + 
                ggplot2::geom_density2d(data = . %>% dplyr::filter(is_symptom==TRUE),
                             alpha=density_alpha,
                             color=symptom_color,
                             adjust=density_adjust)
        }
        
    }
    if(length(root_cells)>0){
        plt <- plt + 
            ggplot2::geom_point(data = . %>% dplyr::filter(is_root==TRUE),
                                color="green",
                                shape=23,
                                alpha=point_alpha*1.5,
                                stroke=1.5,
                                size=3) 
            
    }
    plt <- plt + ggplot2::labs(title = title, subtitle = subtitle)
    #### Reorder layers so that graph and nodes are on top ####
    plt$layers <- c(plt$layers[-graph_layers],plt$layers[graph_layers])
    ## Increase alpha of graph
    plt$layers[[6]]$aes_params$alpha <- .9
    return(
        list(data=cds,
             plot=plt)
    )
}