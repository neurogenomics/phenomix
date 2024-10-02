#' Run pseudotime: subtypes pipeline
#' 
#' Run pseudotime trajectory analysis on a disease and several of its subtypes.
#' @export
#' @importFrom dplyr %>%
#' @examples
#' obj <- get_HPO()
#' id_col <- "id"
#' name_col <- "hpo_name"
run_pseudotime_subtypes <- function(obj,
                                    query,
                                    disease_ids=NULL,
                                    disease_ids_sub=NULL,
                                    title=query,
                                    id_col = "hp_id",
                                    name_col="name",
                                    p2g=HPOExplorer::load_phenotype_to_genes()|>
                                        HPOExplorer::add_disease(),
                                    reduction_ld="umap",
                                    n_genes = 50,
                                    max_disease_ids = 3,
                                    min_gene_sum=0,
                                    future.globals.maxSize = (60*1000)*1024^2,
                                    widths = c(1,.25),
                                    show_plot=TRUE,
                                    force_new=FALSE,
                                    plot_path=tempfile(fileext = ".png"),
                                    data_path=tempfile(fileext = ".csv.gz"),
                                    width = 10, 
                                    height = 16){
    requireNamespace("dplyr")
    disease_name <- NULL;
    force(query)
    
    if(file.exists(plot_path) && 
       file.exists(data_path) &&
       isFALSE(force_new)){
        fig_pseudotime <- magick::image_read(plot_path)
        if(show_plot) methods::show(fig_pseudotime)
        ad_data <- data.table::fread(data_path)
    } else { 
        options(future.globals.maxSize = future.globals.maxSize)
        if(is.null(disease_ids)){
            disease_ids <- c(
                p2g[grepl(query,disease_name,ignore.case = TRUE),]$disease_id|>
                    unique(),
                names(grep(query,obj@meta.data[[name_col]],
                           ignore.case = TRUE, value = TRUE))
            ) |>unique()
        } 
        if(is.null(disease_ids_sub)){
            disease_ids_sub <- disease_ids
        }
        # intersect(colnames(obj), disease_ids)
        # Seurat::NNPlot(o) 
        ## Make disease-level plot
        plot_pseudotime_merged <- plot_pseudotime(
            obj,
            disease_ids = disease_ids,
            # root_cells = disease_ids,
            id_col = id_col,
            title=title,
            merge_trajectories = TRUE,
            p2d = p2g,
            show_plot = FALSE,
            label_roots=TRUE,
            graph_label_size=4,
            cell_size=.5) 
        ## Make phenotype-level subplots
        plot_pseudotime <- plot_pseudotime(
            obj,
            disease_ids = disease_ids_sub,
            max_disease_ids = max_disease_ids,
            id_col = id_col,
            merge_trajectories = FALSE,
            p2d = p2g,
            show_plot = FALSE,
            title_width=20,
            title_col = "disease_name", 
            cell_size=.2)  
        #### Create gene density plot ####
        ## Top AD genes 
        ids <- colnames(obj)[grepl(query,obj@meta.data[[name_col]],
                                   ignore.case = TRUE)]
        genes <- sort(Matrix::rowMeans(Seurat::GetAssayData(obj)[,ids,drop=FALSE]),
                      decreasing = TRUE) 
        if(!is.null(min_gene_sum)) genes <- genes[genes>min_gene_sum]
        plot_feature_density_out <- plot_feature_density(
            obj = obj,
            reduction=reduction_ld,
            joint_only = TRUE,
            method="wkde",
            joint_title = title,
            features = head(names(genes),n_genes)) +
            ggplot2::theme(legend.position = "none")
        ### AD subtype genes
        subtypes <- head(disease_ids,max_disease_ids)
        subtypes_genes <- lapply(stats::setNames(subtypes,
                                                 subtypes), function(x){
            sort(Matrix::rowMeans(Seurat::GetAssayData(obj)[,x,drop=FALSE]),
                 decreasing = TRUE)|> head(1)
        })
        plot_feature_density_out.subtypes <- Nebulosa::plot_density(
            object = obj, 
            reduction = reduction_ld, 
            combine = FALSE,
            features = sapply(subtypes_genes,names)) |>lapply(function(x){
                x+ phenomix::theme_nightlight(legend.position="none")
            }) |> patchwork::wrap_plots() +
            patchwork::plot_layout(axes = "collect", 
                                   axis_titles = "collect",
                                   ncol = 1)
        
        p_pseudo.a <- (
            (plot_pseudotime_merged$plot + 
                 ggplot2::labs(subtitle = paste0(length(disease_ids),
                                                 " disease IDs")) 
            ) |
                (plot_pseudotime$plot +patchwork::plot_layout(ncol = 1))
        ) +patchwork::plot_layout(widths = widths)
        p_pseudo.b <- (
            (
                plot_feature_density_out +
                    ggplot2::labs(title=NULL,
                                  subtitle = paste0("Top ",n_genes,
                                                    " associated genes")) 
            )|
                plot_feature_density_out.subtypes
        ) + patchwork::plot_layout(widths = widths)
        fig_pseudotime <- (p_pseudo.a / p_pseudo.b) + 
            patchwork::plot_annotation(tag_levels = "a") 
        if(show_plot) methods::show(fig_pseudotime)
        #### Save plot ####
        ggplot2::ggsave(plot_path,fig_pseudotime, 
                        width = width, height = height, 
                        dpi = 300)
        #### Save data ####
        ad_data <- subset(plot_pseudotime_merged$plot$data,
                          as.logical(is_root)|as.logical(is_symptom))
        data.table::fwrite(ad_data,data_path)
        gc()
    } 
    return(list(
        data=ad_data,
        plot=fig_pseudotime
    ))
}