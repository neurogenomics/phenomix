#' get_n_descendants
#' 
#' @keywords internal
#' @importFrom ontologyIndex get_descendants
get_n_descendants <- function(ontology,
                              term_ids = NULL) {
    if (is.null(term_ids)) term_ids <- unique(ontology$id)
    n_descends <- lapply(term_ids, function(id) {
        n <- get_ont_level(
            ontology = ontology,
            term_ids = id
        )
        all_descend <- ontologyIndex::get_descendants(
            ontology = ontology,
            roots = id
        )
        n <- length(unique(all_descend))
        return(n)
    }) |>
        `names<-`(term_ids) |>
        unlist()
    return(n_descends)
}
