#' Run enrichment of terms in a DAG ontology
#' 
#' Run enrichment of terms in a DAG ontology.
#' 
#' @export
#' @examples
#' ont <- KGExplorer::get_ontology("hp")
#' obj <- get_HPO()
#' res <- run_dag_enrich(obj, ont, reduction="pca")
run_dag_enrich <- function(obj, 
                           ont,
                           reduction=NULL,
                           cluster_col = NULL,
                           min_hits = 3, 
                           min_offspring = 100,
                           trans_fun=NULL,
                           replace_char=list("."=":",
                                             "_"=":"),
                           q_threshold=.05,
                           top_n=100){
    
    if(is.null(reduction)){
        if(is.null(cluster_col)){
            stopper("Either reduction or cluster_col must be provided.")
        }
        run_dag_enrich_obj(obj=obj, 
                           ont=ont,
                           cluster_col=cluster_col,
                           min_hits=min_hits,
                           min_offspring=min_offspring,
                           replace_char=replace_char,
                           q_threshold=q_threshold)
    } else {
        obsm <- scKirby::get_obsm(obj, keys=reduction, n = 1)
        run_dag_enrich_obsm(obsm=obsm, 
                            ont=ont,
                            min_hits=min_hits,
                            min_offspring=min_offspring,
                            replace_char=replace_char,
                            trans_fun=trans_fun,
                            top_n=top_n,
                            q_threshold=q_threshold)
    }
}