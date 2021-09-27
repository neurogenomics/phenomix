#' Assign ontology groups
#' 
#' Searches for higher-level groups within the ontology
#'  (at depths 2-4 by default), and then assigns each term
#'   in the ontology to one of those groups.
#'  
#' If a given term matches with >1 group at a given depth, 
#' the term will be assigned to the first group arbitrarily.
#' 
#' @param ontology Ontology to search.
#' @param depth_levels Which depth level to assign categories with.
#' @param verbose Print messages.
#' 
#' @return \code{data.table} with group assignments at each depth.
#' 
#' @export
#' @importFrom dplyr %>%
#' @importFrom data.table data.table
get_ontology_groups <- function(ontology,
                                term_ids=NULL,
                                depth_levels=seq(2,4),
                                verbose=TRUE){
    if(is.null(term_ids)) term_ids <- unique(ontology$id)
    depths <- get_n_ancestors(ontology = ontology)
    # meta <- data.table::data.table(id=ontology$id, 
    #                                name=ontology$name)
    # depths <- get_ontology_depths(ontology = ontology)
     # meta$depth <- unname(unlist(depths[meta$id]))
    
    depth_assigns <- lapply(depth_levels, function(d){
        messager("depth =",d,v=verbose)
        top_ids <- names(depths[depths==d])
        top_names <- ontology$name[top_ids]
        top_names <- grep("obsolete|All",top_names,
                          invert = TRUE, value = TRUE)
        all_desc <- lapply(names(top_names), function(root){
            ontologyIndex::get_descendants(ontology = ontology,
                                           roots =  root)
        }) %>% `names<-`(names(top_names))
        
        group_assigns <- lapply(term_ids, function(id){
            group_identity <- lapply(names(all_desc), function(x){
                id %in% all_desc[[x]]
            }) %>% `names<-`(names(all_desc)) %>%
                unlist() 
            group <- top_names[names(group_identity[group_identity][1])]
            return(group)
        }) %>% `names<-`(term_ids) 
        # counts <- table(unname(unlist(group_assigns))) 
        return(group_assigns)
    }) %>% `names<-`(paste0("group_depth",depth_levels)) 
    
    df <- data.frame(do.call("cbind", depth_assigns))
    df <- cbind(id=rownames(df), df)
    #### Remove names within each col (interferes with Seurat) ####
    for(col in colnames(df)){
        df[[col]]  <- unname(unlist(df[[col]]))
    }
    dt <- data.table::data.table(df) 
    return(dt)
}
