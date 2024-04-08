run_dag_enrich_postprocess <- function(RES,
                                       p_threshold,
                                       q_threshold,
                                       sort_by){
    combined_score <- p_value <- z_score <- NULL;
    RES[,p_adjust_all:=stats::p.adjust(p_value, method = "fdr")]
    RES[,combined_score:=(1-p_value)*z_score]
    if(!is.null(p_threshold)) {
        RES <- RES[p_value<p_threshold,]
        if(nrow(RES)==0){
            stopper(paste0("No results left @ p_threshold<",p_threshold))
        }
    }
    if(!is.null(q_threshold)) {
        RES <- RES[p_adjust_all<q_threshold,]
        if(nrow(RES)==0){
            stopper(paste0("No results left @ q_threshold<",q_threshold))
        }
    }
    if(all(names(sort_by) %in% names(RES))){ 
        RES|>
            data.table::setorderv(names(sort_by),
                                  unname(sort_by))
    }
    return(RES)
}