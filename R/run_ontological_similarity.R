#' Run ontologoical similarity
#' 
#' Run Differential ONtological Similairity Analysis (DOSA) on a Seurat object.
#' @inheritDotParams Seurat::FindAllMarkers
#' @export
run_ontological_similarity <- function(obj,
                                       id_col,
                                       ont,
                                       ancestors=unique(ont@elementMetadata$ancestor),
                                       assay_name=paste0(
                                           "ont.",
                                           strsplit(ont@source,",")[[1]][[1]]
                                       ),
                                       group_var="seurat_clusters",
                                       ...){
    force(obj)
    force(id_col)
    force(ont)  
    ontsim <- seurat_to_ontological_similarity(obj=obj,
                                               id_col=id_col,
                                               ont=ont,
                                               ancestors=ancestors,
                                               group_var=group_var)  
    #### Run differential ontological similarity analysis (DOSA) 
    obj[[assay_name]] <- ontsim$assay 
    obj_sim <- obj[,colnames( ontsim$Xsim )]
    messager("Running DOSA on",assay_name)
    Seurat::DefaultAssay(obj_sim) <- assay_name
    markers <- Seurat::FindAllMarkers(obj_sim, 
                                      assay = assay_name,
                                      ...)
    if(nrow(markers)>0){
        top_markers <- markers |> 
            dplyr::filter(p_val_adj < 0.05) |> 
            dplyr::group_by(cluster) |> 
            dplyr::top_n(1, avg_log2FC) |> 
            dplyr::ungroup() |> 
            dplyr::arrange(cluster) |> 
            dplyr::pull(gene) 
    } else {
        top_markers <- markers
    }
    
    #### Return ####
    return(
        list(
            obj_sim=obj_sim,
            ont=ont,
            sim=ontsim$sim,
            Xsim=ontsim$Xsim,
            Xsim_agg=ontsim$Xsim_agg,
            markers=markers,
            top_markers=top_markers
        )
    )
}