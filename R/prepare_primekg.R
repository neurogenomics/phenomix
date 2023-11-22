#' Prepare PrimeKG
#' 
#' Prepare the \href{https://github.com/mims-harvard/PrimeKG}{PrimeKG} dataset
#'  as a \pkg{Seurat} object.
#' @param types Character vector of types of data to return.
#' @param impute Impute missing values in \code{gene_concept} matrix 
#' using the \code{gene_gene} matrix. If \code{TRUE}, the imputed dataset 
#' will be used as a basis of the \code{seurat} object (if selected).
#' @param verbose Print messages.
#' @inheritParams data.table::dcast.data.table
#' @return \code{Seurat} object
#' 
#' @export
prepare_primekg <- function(types=c("gene_concept",
                                    "gene_concept_imputed",
                                    "gene_gene",
                                    "seurat",
                                    "scNLP"),
                            fill = 0,
                            impute=TRUE,
                            verbose = TRUE){
    
    primekg <- data.table::fread("https://dataverse.harvard.edu/api/access/datafile/6180620") 
    primekg[,x_index:=paste0("pkg_",x_index)][,y_index:=paste0("pkg_",y_index)]
    # diseases <- data.table::fread("https://dataverse.harvard.edu/api/access/datafile/6180618")
    # drugs <- data.table::fread("https://dataverse.harvard.edu/api/access/datafile/6180619")
    
    
    res <- list()
    #### Construct gene-concept network ####
    if(any(c("gene_concept","seurat") %in% types)){
        messager("Constructing gene-concept network",v=verbose)
        gc <- primekg[x_type=="gene/protein" & y_type!="gene/protein"][,dummy:=1]
        res[["gene_concept"]] <- data.table::dcast.data.table(
            gc,
            formula = x_index ~ y_index,
            value.var = "dummy",
            fill = fill,
            fun.aggregate = "mean",
            na.rm = TRUE
        ) |> scKirby::to_sparse() 
        remove(gc)
    }
   
    #### Construct gene-gene network ####
    if("gene_gene" %in% types){
        messager("Constructing gene-gene network",v=verbose)
        gg <- primekg[x_type=="gene/protein" & y_type=="gene/protein"][,dummy:=1]
        res[["gene_gene"]] <- data.table::dcast.data.table(
            gg,
            formula = x_index ~ y_index,
            value.var = "dummy",
            fill = fill,
            fun.aggregate = "mean",
            na.rm = TRUE
        ) |> scKirby::to_sparse() 
        remove(gg)
    }
    
    #### Impute ####
    if(isTRUE(impute)){
        Xnet = cbind(O=rep(1,nrow(res[["gene_gene"]])),
                     res[["gene_gene"]])
        Xi <- run_imputation(X = res[["gene_concept"]], 
                             Xnet = Xnet,
                             cores = -1)
        res[["gene_concept_imputed"]] <- Xi$Network |> scKirby::to_sparse()
        remove(Xi)
    }
    
    #### Construct Seurat object ####
    if("seurat" %in% types){
        obs <- gc[,grep("^y_",names(gc), value = TRUE), with = FALSE] |> unique()
        obj <- list(data = list(imputed=res$gene_concept_imputed,
                                raw=res$gene_concept),
                    obs = data.frame(obs, row.names = obs$y_index))
        res[["seurat"]] <- scKirby::process_seurat(obj = obj,
                                                   nfeatures = NULL) 
        # Seurat::DimPlot(res[["seurat"]],group.by = "y_type")
        if("scNLP" %in% types){
            res[["scNLP"]] <- scNLP::plot_tfidf(
                obj = res[["seurat"]],
                point_palette = pals::kovesi.cyclic_mygbm_30_95_c78_s25(length(unique(res[["seurat"]]$seurat_clusters))),
                point_alpha = .1,
                terms_per_cluster = 2,
                size = "nFeature_imputed",
                reduction = "umap",
                label_var = "y_name")
        }
    }
    #### Return ####
    return(res)
}
