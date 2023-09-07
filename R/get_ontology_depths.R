get_ontology_depths <- function(ontology,
                                term_ids = NULL,
                                as_depths = TRUE,
                                workers = NULL) {
    if (is.null(term_ids)) term_ids <- unique(ontology$id)
     
    # lvls <- parallel::mclapply(term_ids, function(id){
    #     ancests <- ontologyIndex::get_ancestors(ontology = ontology,
    #                                             terms = id)
    #     lvls <- get_ontology_levels(ontology = ontology,
    #                                 term_ids = ancests)
    #     lvl <- min(lvls,na.rm = TRUE)+1
    #     return(lvl)
    # }, mc.cores = cores$worker_cores) |>
    #     `names<-`(term_ids) |> unlist()
    lvls <- get_ontology_levels(
        ontology = ontology,
        term_ids = term_ids
    )

    if (as_depths) {
        # Invert so that larger numbers means deeper into the ontology
        # (ie more granular terms)
        depths <- max(lvls, na.rm = TRUE) - lvls
        return(depths)
    } else {
        return(lvls)
    }
}
