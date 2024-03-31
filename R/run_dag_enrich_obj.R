run_dag_enrich_obj <- function(obj, 
                               ont,
                               cluster_col,
                               min_hits,
                               min_offspring,
                               replace_char,
                               q_threshold=.05){
    z_score <- NULL;
    messager("Running dag_enrich_on_offsprings.")
    clusters <- sort(unique(obj[[cluster_col]][,1]))
    lapply(seq(length(clusters)), function(i){
        nms <- colnames(obj)[obj[[cluster_col]][,1]==i]
        nms <- replace_char_fun(nms,replace_char)
        nms <- intersect(nms,ont@terms)
        if(length(nms)<min_hits) {
            messager("Not enough terms found for cluster:",i)
            return(NULL)
        }
        res <- simona::dag_enrich_on_offsprings(dag = ont, 
                                                terms = nms,
                                                min_hits = min_hits,
                                                min_offspring = min_offspring) |>
            data.table::data.table(key="term")
        res <- cbind(group=clusters[i],res)
        res[,combined_score:=(1-p_adjust)*z_score]
        if(is.null(q_threshold)) return(res)
        res_sig <- res[p_adjust<q_threshold,][!is.na(combined_score)]|>
            data.table::setorderv(c("combined_score"),c(-1))
        res_sig
    }) |> data.table::rbindlist()
}