#' Seurat to ontological similarity
#' 
#' Takes a Seurat object and computes the ontological similarity of 
#' each sample (column) given some ontology that it can be mapped onto.
#' @inheritDotParams simona::term_sim
#' @export
seurat_to_ontological_similarity  <- function(obj,
                                              id_col,
                                              ont,
                                              ancestors,
                                              group_var=NULL,
                                              as_sparse=TRUE,
                                              return_assay=TRUE,
                                              ...){
    if(!id_col %in% colnames(obj@meta.data)){
        stop("id_col not found in obj@metadata.")
    }
    ## Compute similarity of each term to each ancestor term in HPO 
    matched_meta <- obj@meta.data[
        simona::dag_has_terms(ont,terms = obj@meta.data[[id_col]]) |
        simona::dag_has_terms(ont,terms = colnames(obj)),
    ]
    sim <- simona::term_sim(ont, 
                            terms = unique(c(ancestors,
                                             matched_meta[[id_col]]),
                            ...
                            )
    ) 
    # for(x in unique(ancestors)){
    #   message("Adding ontological similarity for term: ",x)
    #   obj@meta.data[[paste0("sim.",x)]] <- sapply(obj$hp_id, function(x){
    #   if(x %in% colnames(sim)){
    #     sim[x,ancestors[1]]
    #   } else {
    #     NA
    #   }})
    # }
    ### Map names back onto original IDs in obj
    og_ids <- stats::setNames(matched_meta[[id_col]],
                              rownames(matched_meta)
    )[matched_meta[[id_col]] %in% colnames(sim)]
    Xsim <- sim[ancestors,unname(og_ids)]
    colnames(Xsim) <- names(og_ids)
    ## Add as a new assay
    if(isTRUE(return_assay)){
        extra_ids <- setdiff(colnames(obj),colnames(Xsim))
        if(length(extra_ids)>0){
            Xextra <- matrix(NA,
                             nrow = nrow(Xsim),
                             ncol = length(extra_ids)
            )|>`colnames<-`(extra_ids)
            Xsim <- cbind(Xsim,Xextra)
        }
        if(as_sparse) Xsim <- methods::as(Xsim,"sparseMatrix")
        ont_assay <- Seurat::CreateAssayObject(data = Xsim)
    } else{
        ont_assay <- NULL
    } 
    #### Compute mean similarity to each ancestor by cluster
    if(!is.null(group_var)){
        Xsim_agg<- orthogene:::aggregate_rows(
            X = Matrix::t(Xsim),
            agg_fun="mean",
            groupings = unname(obj@meta.data[colnames(Xsim),][[group_var]]),
            as_DelayedArray=FALSE)
        # heatmap(as.matrix(Xsim_agg))
    } else {
        Xsim_agg <- NULL
    } 
    return(list(
        sim=sim,
        Xsim=Xsim,
        Xsim_agg=Xsim_agg,
        assay=ont_assay
    ))
}