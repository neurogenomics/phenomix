run_dag_enrich_obsm <- function(obj, 
                                ont,
                                id_col,
                                reduction,
                                min_hits,
                                min_offspring,
                                replace_char,
                                top_n,
                                trans_fun,
                                value_threshold,
                                p_threshold,
                                q_threshold,
                                sort_by,
                                on_offspring=TRUE){
    
    obsm <- scKirby::get_obsm(obj, keys=reduction, n = 1)
    # ont_prefixes <- unique(stringr::str_split(ont@terms,pattern = ":",
    #                                           simplify = TRUE)[,1])
    # obsm <- obsm[grepl(paste(paste0("^",ont_prefixes),collapse = "|"),
    #                    rownames(obsm)),]
    if(!is.null(trans_fun)) obsm <- trans_fun(obsm)
    if(nrow(obsm)==0) stopper("No matching terms found in obsm.")
    id_dict <- if(is.null(id_col)) {
        stats::setNames(colnames(obj),
                        colnames(obj))
    } else {
        if(!id_col %in% colnames(obj@meta.data)){
            messager("Could not find",paste0("id_col=",shQuote(id_col)),
                     "Defaulting to colnames(obj).")
            stats::setNames(colnames(obj),
                            colnames(obj))
        }
        stats::setNames(obj[[id_col]][,1],
                        rownames(obj[[id_col]]))
    }
    messager("Running dag_enrich_on_offsprings.")
    BPPARAM <- KGExplorer::set_cores()
    RES <- BiocParallel::bplapply(stats::setNames(seq(ncol(obsm)),
                                                  colnames(obsm)), 
                                  BPPARAM = BPPARAM, 
                                  function(i){
        val <- obsm[,i]
        if(!is.null(value_threshold)) val <- val[val>value_threshold]
        nms <- id_dict[names(tail(sort(val),top_n))]|>unname()
        nms <- map_id_sep(nms,replace_char)
        nms <- intersect(nms,ont@terms)
        if(length(nms)<min_hits) {
            messager("Not enough terms found for",colnames(obsm)[i],
                     v=!BPPARAM$progressbar)
            return(NULL)
        }
        if(on_offspring){
            res <- simona::dag_enrich_on_offsprings(dag = ont, 
                                                    terms = nms,
                                                    min_hits = min_hits,
                                                    min_offspring = min_offspring)
        } else {
            res <- simona::dag_enrich_on_items(dag = ont, 
                                               items = nms,
                                               min_hits = min_hits,
                                               min_items = min_offspring)
        }
        res <- data.table::data.table(res, key="term")
        res <- cbind(group=colnames(obsm)[i],res)
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
