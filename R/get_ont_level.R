get_ont_level = function(ontology,
                         term_ids) {
    # Written by Bobby Gordon-Smith
    children = unique(setdiff(unlist(ontology$children[term_ids]), term_ids))
    if (length(children) == 0) {
        return(0)
    } else {
        return(1 + get_ont_level(ontology,children)) #<- recursion..
    }
}