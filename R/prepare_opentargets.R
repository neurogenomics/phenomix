#' Prepare OpenTargets
#'
#' Prepare the OpenTargets gene-by-disease associations as a 
#' \link{Seurat} object.
#' @inheritParams scKirby::process_seurat
#' @inheritDotParams scKirby::process_seurat
#' @import KGExplorer
#' @export
prepare_opentargets <- function(data_type="associationByOverallDirect",
                                formula = "approvedSymbol ~ diseaseId",
                                nfeatures=NULL,
                                default_assay = "score",
                                vars.to.regress = paste0("nFeature_",
                                                         default_assay),
                                save_path=NULL,
                                force_new=FALSE,
                                ...){
    targetId <- id <- NULL;
    
    if(file.exists(save_path) && !force_new){
        messager("Loading precomputed data:",save_path)
        return(readRDS(save_path))
    } 
    #### Get gene-disease relationship data ####
    d <- KGExplorer::get_opentargets(data_type=data_type)
    #### Get sample metadata ####
    obs <- KGExplorer::get_opentargets(data_type = "diseases")
    obs$db <- stringr::str_split(obs$id,"_",simplify = TRUE)[,1]
    messager("Mapping xref IDs.")
    for(x in unique(obs$db)){
        suppressWarnings(
            map_xref(obs,prefix=x,verbose=FALSE)    
        )
    }
    #### Get feature metadata ####
    var <- KGExplorer::get_opentargets(data_type = "targets")
    var <- var[id %in% unique(d$targetId)] |>
        data.table::setkeyv("id")
    #### Convert ensembl IDs to gene symbols ####
    d[,approvedSymbol:=var[targetId]$approvedSymbol]
    #### Convert to sparse matrix ####
    messager("Constructing matrix:",formula)
    X <- data.table::dcast.data.table(d,
                                      formula = as.formula(formula),
                                      value.var = "score",
                                      fun.aggregate = mean, 
                                      fill = 0,
                                      na.rm = TRUE) |>
        scKirby::to_sparse()
    #### Subset by intersecting IDs ####
    ids <- intersect(obs$id, d$diseaseId)
    obs <- obs[id %in% ids,]
    X <- X[,ids]
    var <- var[approvedSymbol %in% rownames(X),.SD[1],by="approvedSymbol"]
    #### Remove columns in obs/var that are lists ####
    obs <- obs[,names(obs)[!sapply(obs,is.list)], with=FALSE]
    #### Convert list columns to characters ####
    obs <- KGExplorer::unlist_dt(obs, drop=TRUE)
    var <- KGExplorer::unlist_dt(var, drop=TRUE)
    #### Convert to Seurat ####
    messager("Processing as Seurat object.")
    obj <- scKirby::process_seurat(
        obj = list(data = list(X)|> `names<-`(default_assay), 
                   obs = data.frame(obs,row.names = obs$id),
                   var = data.frame(var,row.names = var$approvedSymbol)),
        nfeatures = nfeatures,
        vars.to.regress = vars.to.regress,
        default_assay = default_assay,
        ...)
    #### Save ####
    if(!is.null(save_path)) cache_save(obj,save_path)
    #### Return ####
    return(obj)
}