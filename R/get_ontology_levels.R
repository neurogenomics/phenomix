get_ontology_levels <- function(ontology,
                                term_ids = NULL) {
    if (is.null(term_ids)) term_ids <- unique(ontology$id)
    lvls <- lapply(term_ids, function(id) {
        get_ont_level(
            ontology = ontology,
            term_ids = id
        )
    }) |>
        `names<-`(term_ids) |>
        unlist()
    return(lvls)
}
