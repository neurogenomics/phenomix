#' Map phenotypes 
#' 
#' Map phenotypes onto a standardised ontology.
#' @param phenotypes A character vector of phenotype terms.
#' @param ont An ontology of class \code{ontology_index}.
#' @returns Named vector mapping the phenotypes onto an ontology.
#' 
#' @export
#' @import data.table
#' @examples
#' meta <- MungeSumstats::find_sumstats()
#' phenotypes <- unique(meta$trait)
#' pmap <- map_phenotypes(phenotypes=phenotypes) 
map_phenotypes <- function(phenotypes,
                           ont = KGExplorer::get_ontology("hp")){ 
    requireNamespace("HPOExplorer")
    
    pfilt <- phenotypes[tolower(phenotypes) %in% tolower(ont$name)] 
    ont$name2 <- stats::setNames(names(ont$name), tolower(ont$name)) 
    pmap <- data.table::data.table(
        phenotype=pfilt,
        ont_id=ont$name2[tolower(pfilt)],
        ont_name=ont$name[ont$name2[tolower(pfilt)]],
        key = "phenotype")
    return(pmap)
}
