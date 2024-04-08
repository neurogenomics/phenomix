#' Plot top cell types
#' 
#' Plot the top cell types per cluster.
#' @export
#' @examples
#' obj <- get_HPO()
#' out <- plot_top_celltypes(obj)
plot_top_celltypes <- function(obj,
                               types=c("radial","pie","reduction"),
                               reduction="umap",
                               prefix="^q[.]",
                               cluster_vars=c("seurat_clusters"),
                               color_col=cluster_vars[1],
                               label_col=cluster_vars[1],
                               enriched_proportion_threshold=.2,
                               label_alpha=.9,
                               label_size=3,
                               label_color="white"){
    requireNamespace("ggplot2")
    variable <- NULL;
    
    obj <- add_cluster_colors(obj)
    cluster_colors <- add_cluster_colors(obj, return_dict = TRUE)
    #### Assign cell type ancestors ####
    # cl <- KGExplorer::get_ontology("cl")
    # cl <- KGExplorer::filter_ontology(cl,
    #                                   remove_terms = c("continuant",
    #                                                    "electrically active cell",
    #                                                    "nucleate cell",
    #                                                    "stuff accumulating cell",
    #                                                    "precursor cell"),
    #                                   keep_descendants = c("cell"))
    # cl <- KGExplorer::add_ancestors(cl,lvl = 1, force_new = TRUE, i = 2)
    # cl_dt <- data.table::data.table(cl@elementMetadata, key="name")
    # meta_mean[,top_ancestor:=cl_dt[top_celltype]$ancestor_name[2], by=c(cluster_vars)]
    # table(meta_mean$top_ancestor, useNA = "always")
    
    #### Get cell types ####
    celltypes <- grep(prefix,colnames(obj@meta.data), value = TRUE)
    if(length(celltypes)==0)stopper(paste0(
        "No cell type columns found with prefix=",
        shQuote(prefix)
    ))
    {
        meta <- data.table::data.table(obj@meta.data, keep.rownames = "rn")
        cluster_vars <- intersect(cluster_vars, names(meta))
        embed <- scKirby::get_obsm(obj,
                                   keys = reduction, 
                                   n=1)|>
            data.table::data.table(keep.rownames = "rn") 
        embed_cols <- names(embed)[-1]
        for(e in embed_cols){
            if(e %in% names(meta)) meta <- meta[,-c(e), with=FALSE]
        }
        meta <- merge(meta,embed, by="rn")
    }
    meta_melt <- (
        meta|>
            data.table::melt.data.table(measure.vars = celltypes,
                                        value.name = "q")
    )[!is.na(q),]
    meta_melt[,celltype:=gsub("[.]"," ",
                              gsub(prefix,"",variable))][is.na(q),q:=1] |>
        data.table::setorderv("q")
    #### Aggregate at cluster level ####
    meta_melt[,c("total_sig_count",
                 "total_sig_celltypes",
                 "total_count",
                 embed_cols):=list(sum(q<0.05),
                                 data.table::uniqueN(celltype[q<0.05]),
                                 .N,
                                 mean(get(embed_cols[1])),
                                 mean(get(embed_cols[2]))
                                 ),
              by=c(cluster_vars)] 
    meta_agg <- meta_melt[,list(
        celltype_sig_count=sum(q<0.05),
        celltype_count=.N,
        total_sig_count=unique(total_sig_count),
        total_count=unique(total_count),
        mean_q=mean(q),
        umap_1=mean(get(embed_cols[1])),
        umap_2=mean(get(embed_cols[2]))
    ),
    by=c(unique(c(cluster_vars,label_col,"celltype")))]|>
        data.table::setorderv(c(cluster_vars,"celltype_sig_count","mean_q"),
                              c(rep(1,length(cluster_vars)),-1,1))
    
    out <- list()
    if("radial" %in% types){
        out[["radial"]][["data"]] <- meta_agg
        out[["radial"]][["plot"]] <- (
            ggplot2::ggplot(
                data=meta_agg, 
                ggplot2::aes(x=seurat_clusters,  
                             y=celltype_sig_count, 
                             fill=celltype)) +
                ggplot2::geom_bar(stat="identity", position = "fill") +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                                   hjust = 1)) +
                ggplot2::labs(title = "Cell type associations",x=NULL,y=NULL) +
                ggplot2::theme(legend.position = "none") +
                ggplot2::theme_minimal() +
                ggplot2::coord_polar("x", start=0) +
                theme_nightlight()
        ) |>
            shrink_legend(ncol=4)
    }
    
    if("pie" %in% types){
        out[["pie"]][["data"]] <- meta_agg
        out[["pie"]][["plot"]] <- (
            ggplot2::ggplot(meta_agg, 
                            ggplot2::aes(x=seurat_clusters, 
                                         y=celltype_sig_count, 
                                         fill=celltype)) +
                ggplot2::geom_bar(stat="identity") +
                ggplot2::facet_wrap(facets = as.formula(
                    paste("~",
                          paste(unique(c(color_col,label_col)),collapse=" + ")
                    )
                ), 
                scales = "free",
                nrow=6) + 
                ggplot2::coord_polar("y", start=0) +
                theme_nightlight() +
                ggplot2::labs(title="seurat_clusters",x=NULL,y=NULL) 
        )|>
            shrink_legend()
            
    }
    if("reduction" %in% types){
        meta_agg[,enriched_proportion:=(celltype_sig_count/celltype_count),
                 by=c(cluster_vars,"celltype")]|>
            data.table::setorderv("enriched_proportion",-1)
        plot_dat <- meta_agg[!is.na(get(label_col)),.SD[1], by=c(cluster_vars)] 
        
        out[["reduction"]][["data"]] <- plot_dat
        out[["reduction"]][["plot"]] <- 
            Seurat::DimPlot(obj, 
                            group.by = color_col) +
            Seurat::NoLegend() + 
            ggplot2::scale_color_manual(values = cluster_colors) +
            ggrepel::geom_label_repel(
                data=plot_dat[enriched_proportion>enriched_proportion_threshold],
                ggplot2::aes(x=!!ggplot2::sym(embed_cols[1]), 
                             y=!!ggplot2::sym(embed_cols[2]),
                             label=paste0(!!ggplot2::sym(label_col),
                                          "\n[",celltype,"]=",
                                          round(enriched_proportion,2)),
                             fill=!!ggplot2::sym(color_col),
                             size=enriched_proportion
                ), 
                alpha=label_alpha,
                size=label_size,
                color=label_color,
                arrow = ggplot2::arrow(type="closed",
                                       length = ggplot2::unit(0.1, "inches")),
                force = 2,
                box.padding = 2,
                min.segment.length = 0,
                max.overlaps=100) +
            ggplot2::scale_fill_manual(values = cluster_colors) +
            theme_nightlight()
    }
    #### Return ####
    return(out)
}