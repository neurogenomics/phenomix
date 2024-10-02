#' Plot ontological similarity
#' 
#' Plot ontological similarity of cells to a set of terms in an ontology.
#' @inheritDotParams Seurat::FeaturePlot
#' @export
plot_ontological_similarity <- function(obj,
                                        id_col,
                                        ont,
                                        ### Use pre-computed DOSA results 
                                        run_ontological_similarity_out=NULL,
                                        ### Args for running new DOSA
                                        ancestors=unique(ont@elementMetadata$ancestor),
                                        group_var="seurat_clusters",
                                        top_n=9,
                                        min.cutoff = "q90",
                                        show_plot=TRUE,
                                        ...){
    if(!is.null(run_ontological_similarity_out)){ 
        messager("Using pre-computed DOSA results.")
        out <- run_ontological_similarity_out
    } else {
        out <- run_ontological_similarity(
            obj=obj,
            id_col=id_col,
            ont=ont,
            ancestors=ancestors,
            group_var=group_var) 
    }
    #### Plot most differentially similar ancestors ####
    terms_to_use <- if(length(out$top_markers)==0){
        messager("No term identified from DOSA. Using ancestors instead.")
         intersect(ancestors,colnames(out$sim))
    } else {
        unique(out$top_markers)
    }
    term_map <- KGExplorer::map_ontology_terms(ont = ont,
                                               keep_order = FALSE,
                                               terms = terms_to_use)|>
        utils::head(top_n) 
    term_map <- term_map[!is.na(names(term_map))]
    ### Order ancestors by their clustered similarity
    d <- stats::dist(out$sim[names(term_map),])
    hc <- stats::hclust(d) 
    term_map <- term_map[hc$labels[hc$order]]
    
    fo <- Seurat::FeaturePlot(out$obj_sim,
                              order = TRUE,  
                              min.cutoff = min.cutoff,
                              features = names(term_map),
                              ...) 
    fig_ontmarkers <- lapply(seq(length(term_map)), function(i){
        fo_i <- fo[[i]] +
            ggplot2::labs(title=NULL,
                          subtitle=paste(term_map[[i]],
                                         names(term_map)[[i]],sep="\n"),
                          color="Ontological\nsimilarity"
            )+
            phenomix::theme_nightlight() +
            phenomix::scale_color_nightlight()
        if(i>1){
            fo_i <- fo_i + ggplot2::theme(legend.position = "none")
        }
        fo_i
    })|>patchwork::wrap_plots() +
        patchwork::plot_layout(axes = "collect",
                               axis_titles = "collect",
                               guides = "collect"
        )
    
    if(isTRUE(show_plot)) methods::show(fig_ontmarkers)
    #### Return ####
    return(
        list(plot=fig_ontmarkers,
             run_ontological_similarity=out,
             hclust=hc,
             term_map=term_map
        )
    )
}