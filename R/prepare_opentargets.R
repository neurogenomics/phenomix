#' Prepare OpenTargets
#'
#' Prepare the OpenTargets gene-by-disease associations as a 
#' \link{Seurat} object.
#' @inheritParams scKirby::process_seurat
#' @inheritDotParams scKirby::process_seurat
#' @import KGExplorer
#' @export
prepare_opentargets <- function(data_type="associationByDatasourceIndirect",#"associationByOverallDirect",
                                formula = "approvedSymbol ~ id",
                                nfeatures=NULL,
                                default_assay = "score",
                                vars.to.regress = NULL, #paste0("nFeature_",default_assay),
                                fill=0,
                                save_path=NULL,
                                force_new=FALSE,
                                ...){
    targetId <- id <- NULL;
    
    if(!is.null(save_path) && 
       file.exists(save_path) && 
       isFALSE(force_new)){
        messager("Loading precomputed data:",save_path)
        return(readRDS(save_path))
    } 
    #### Get gene-disease relationship data ####
    d <- KGExplorer::get_opentargets(data_type = data_type,
                                     force_new = force_new)
    #### Get sample metadata ####
    obs <- KGExplorer::get_opentargets(data_type = "diseases",
                                       force_new = force_new)
    obs$db <- stringr::str_split(obs$id,"_",simplify = TRUE)[,1]
    obs|> data.table::setnames("description","definition")
    messager("Mapping xref IDs.")
    for(x in unique(obs$db)){
        suppressWarnings(
            map_xref(dat = obs,
                     prefix=x,
                     verbose=FALSE)    
        )
    }
    #### Get feature metadata ####
    ## Mapping with orthogene zero difference when converting to gene symbols
    # if(isTRUE(run_map_genes)){
    #     var <- KGExplorer::get_opentargets(data_type = "targets")
    #     var <- var[id %in% unique(d$targetId)] 
    #     gene_map <- orthogene::map_genes(genes = unique(var$approvedSymbol))
    #     gene_map <- data.table::data.table(gene_map)[,.SD[1], by="input"]
    #     gene_map[,approvedSymbol:=data.table::fcoalesce(name,input)]
    #     data.table::setkeyv(gene_map,"input")
    #     var[,approvedSymbol_standard:=data.table::fcoalesce(
    #         gene_map[approvedSymbol]$name,approvedSymbol)] 
    #     data.table::setkeyv(var,"id")
    #     #### Convert ensembl IDs to gene symbols ####
    #     d[,approvedSymbol:=var[targetId]$approvedSymbol_standard]   
    # } else{
        #### Get feature metadata ####
        var <- KGExplorer::get_opentargets(data_type = "targets", 
                                           force_new = force_new)
        var <- var[id %in% unique(d$targetId)] |>
            data.table::setkeyv("id")
        #### Convert ensembl IDs to gene symbols ####
        d[,approvedSymbol:=var[targetId]$approvedSymbol]
    # }
    #### Add IDs ####
    if("datasourceId" %in% names(d)){
        d[,id:=paste(map_id_sep(diseaseId),datasourceId,sep=".")]
    } else {
        d[,id:=map_id_sep(diseaseId)]
    }
    #### Convert to sparse matrix ####
    messager("Constructing matrix:",formula)
    ### Chunk by file to avoid memory issues
    chunks <- unique(d$file)
    BPPARAM <- KGExplorer::set_cores()
    X <- BiocParallel::bplapply(chunks,
                                 BPPARAM = BPPARAM,
                                 function(f){
        data.table::dcast.data.table(d[file==f,],
                                     formula = as.formula(formula),
                                     value.var = "score",
                                     fun.aggregate = mean, 
                                     fill = fill,
                                     na.rm = TRUE) |>
            KGExplorer::dt_to_matrix(as_sparse = TRUE)
    })
    messager("Merging",length(X),"sparse matrices.")
    X <- Reduce(Seurat::RowMergeSparseMatrices,X)
    #### Report matrix stats ####
    messager("Merged matrix:",nrow(X),"x",ncol(X))
    ## Check gene sparsity 
    ngenes <- Matrix::colSums(X>0)
    messager(sum(ngenes<2),"/",length(ngenes),
             paste0("(",round(sum(ngenes<2)/length(ngenes)*100,1),"%)"),
             "traits have <2 genes.")
    #### Subset by intersecting IDs ####
    obs <- merge(
        unique(d[,-c("targetId","score","evidenceCount","approvedSymbol")]),
          obs[,-c("file")],
          by.x="diseaseId",
         by.y="id",
        all.x=TRUE)
    X <- X[,obs$id]
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
    KGExplorer::cache_save(obj,save_path)
    #### Return ####
    return(obj)
}