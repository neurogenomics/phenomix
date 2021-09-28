#' get_n_ancestors
#' 
#' @keywords internal
#' @importFrom ontologyIndex get_ancestors
get_n_ancestors <- function(ontology,
                            term_ids = NULL) {
    if (is.null(term_ids)) term_ids <- unique(ontology$id)
    n_ancestors <- lapply(term_ids, function(x) {
        all_terms <- ontologyIndex::get_ancestors(
            ontology = ontology,
            terms = x
        )
        return(length(all_terms))
    }) %>%
        `names<-`(term_ids) %>%
        unlist()
    return(n_ancestors)
}
