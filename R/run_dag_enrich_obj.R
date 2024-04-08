run_dag_enrich_obj <- function(obj, 
                               ont,
                               id_col,
                               cluster_col,
                               min_hits,
                               min_offspring,
                               replace_char,
                               p_threshold,
                               q_threshold,
                               sort_by){
    messager("Running dag_enrich_on_offsprings.")
    clusters <- sort(unique(obj[[cluster_col]][,1]))
    all_ids <- if(is.null(id_col)) {
        colnames(obj)
    } else {
        if(!id_col %in% colnames(obj@meta.data)){
            messager("Could not find",paste0("id_col=",shQuote(id_col)),
                     "Defaulting to colnames(obj).")
            colnames(obj)
        }
        obj@meta.data[[id_col]]
    }
    BPPARAM <- KGExplorer::set_cores()
    RES <- BiocParallel::bplapply(seq(length(clusters)),
                           BPPARAM = BPPARAM, 
                           function(i){ 
        nms <- all_ids[obj@meta.data[[cluster_col]]==i]
        nms <- map_id_sep(nms,replace_char)
        nms <- intersect(nms,ont@terms)
        if(length(nms)<min_hits) {
            warning("Not enough terms found for cluster:",i)
            return(NULL)
        }
        res <- simona::dag_enrich_on_offsprings(dag = ont, 
                                                terms = nms,
                                                min_hits = min_hits,
                                                min_offspring = min_offspring) |>
            data.table::data.table(key="term")
        res <- cbind(group=clusters[i],res)
        res$input_ids <- paste(nms,collapse = ";")
        res$input_names <- paste(KGExplorer::map_ontology_terms(ont = ont, 
                                                                terms = nms,
                                                                to = "name",
                                                                verbose = FALSE),
                                 collapse = ";")
        return(res)
    }) |> data.table::rbindlist(fill=TRUE) 
    run_dag_enrich_postprocess(RES=RES,
                               p_threshold=p_threshold,
                               q_threshold=q_threshold,
                               sort_by=sort_by)
}
