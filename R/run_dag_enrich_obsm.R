run_dag_enrich_obsm <- function(obsm, 
                                ont,
                                min_hits,
                                min_offspring,
                                replace_char,
                                top_n=100,
                                trans_fun=NULL,
                                q_threshold=.05){
    z_score <- NULL;
    
    ont_prefixes <- unique(stringr::str_split(ont@terms,pattern = ":",
                                              simplify = TRUE)[,1])
    obsm <- obsm[grepl(paste(paste0("^",ont_prefixes),collapse = "|"),
                       rownames(obsm)),]
    if(!is.null(trans_fun)) obsm <- trans_fun(obsm)
    if(nrow(obsm)==0) stopper("No matching terms found in obsm.")
    messager("Running dag_enrich_on_offsprings.") 
    lapply(stats::setNames(seq(ncol(obsm)),
                           colnames(obsm)), 
           function(i){
        val <- obsm[,i]
        nms <- names(tail(sort(val),top_n))
        nms <- replace_char_fun(nms,replace_char)
        nms <- intersect(nms,ont@terms)
        if(length(nms)<min_hits) {
            messager("Not enough terms found for",colnames(obsm)[i])
            return(NULL)
        }
        res <- simona::dag_enrich_on_offsprings(dag = ont, 
                                                terms = nms,
                                                min_hits = min_hits,
                                                min_offspring = min_offspring) |>
            data.table::data.table(key="term")
        res <- cbind(group=colnames(obsm)[i],res)
        res[,combined_score:=(1-p_adjust)*z_score]
        if(is.null(q_threshold)) return(res)
        res_sig <- res[p_adjust<q_threshold,][!is.na(combined_score)]|>
            data.table::setorderv(c("combined_score"),c(-1))
        res_sig
    }) |> data.table::rbindlist()
}