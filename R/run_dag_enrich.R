#' Run enrichment of terms in a DAG ontology
#' 
#' Run enrichment of terms in a DAG ontology.
#' 
#' @param value_threshold Minimum weight to include a term 
#' (after applying \code{truns_fun}).
#'  Only used when \code{reduction} is provided.
#' @param on_offspring Use the function \link[simona]{dag_enrich_on_offsprings}
#'  (TRUE) or \link[simona]{dag_enrich_on_items} (FALSE).
#' @export
#' @examples
#' ont <- KGExplorer::get_ontology("hp")
#' obj <- get_HPO()
#' res <- run_dag_enrich(obj, ont, reduction="pca")
run_dag_enrich <- function(obj, 
                           ont,
                           id_col=NULL,
                           reduction=NULL,
                           cluster_col = NULL,
                           min_hits = 3, 
                           min_offspring = 100,
                           min_depth=NULL,
                           trans_fun=NULL,
                           replace_char=list("."=":",
                                             "_"=":"),
                           value_threshold=NULL,
                           p_threshold=.05,
                           q_threshold=.05,
                           top_n=100,
                           sort_by= c("log2_fold_enrichment"=-1),
                           on_offspring=TRUE
                           ){
    if(!is.null(min_depth)){
        KGExplorer::filter_ontology()
    }
    
    if(is.null(reduction)){
        if(is.null(cluster_col)){
            stopper("Either reduction or cluster_col must be provided.")
        }
        run_dag_enrich_obj(obj=obj, 
                           ont=ont,
                           id_col=id_col,
                           cluster_col=cluster_col,
                           min_hits=min_hits,
                           min_offspring=min_offspring,
                           replace_char=replace_char,
                           p_threshold=p_threshold,
                           q_threshold=q_threshold,
                           sort_by=sort_by,
                           on_offspring=on_offspring)
    } else {
        run_dag_enrich_obsm(obj=obj, 
                            ont=ont,
                            id_col=id_col,
                            reduction=reduction,
                            min_hits=min_hits,
                            min_offspring=min_offspring,
                            replace_char=replace_char,
                            trans_fun=trans_fun,
                            top_n=top_n,
                            value_threshold=value_threshold,
                            p_threshold=p_threshold,
                            q_threshold=q_threshold,
                            sort_by=sort_by,
                            on_offspring=on_offspring)
    }
}