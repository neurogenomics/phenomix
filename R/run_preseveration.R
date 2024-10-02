#' Run preservation
#'
#' Assess the preservation of high-dimensional structure in 
#' mid-dimensional and low-dimensional space.
#' @param obj Seurat object.
#' @param reduction List of reductions to use for preservation.
#' @param id_types List of id types to use for preservation.
#' @param k.param Number of nearest neighbours to use for preservation.
#' @param dims List of dimensions to use for preservation.
#' @param cor_fun Function to use for correlation.
#' @param max_ids Maximum number of ids to use for preservation.
#' @return data.table with preservation results.
#' @export
#' @examples
#' obj <- get_HPO()
#' obj$id[1:100] <- obj$id[1] 
#' res <- run_preservation(obj, id_types="id")
run_preservation <- function(obj,
                             reduction=list(md="pca",
                                            ld="umap"),
                             id_types = c("hp_id","mondo_id","efo_id"),
                             k.param =30,
                             dims=list(md=1:10,
                                       ld=1:2),
                             cor_fun=run_cor,
                             max_ids=10000,
                             save_path=NULL,
                             force_new=FALSE
                             ){
    id_type <- id <- n <- diff_cor_hd <- diff_cor_md <- diff_cor_ld <- NULL;
    if(!is.null(save_path) &&
       file.exists(save_path) && 
       isFALSE(force_new)){
        messager("Loading preservation results from:",save_path)
        return(readRDS(save_path))
    } else{
        id_types <- intersect(id_types, colnames(obj@meta.data))
        if(length(id_types)==0) stopper("No id_types found in obj metadata.")
        ## Seurat returns KNN graphs in a very weird way.
        ## Contrary to convention, larger values mean smaller distances (greater overlap).
        ## Shared NN is normalised from 0-1, where 1 indicates perfect overlap.
        if(!"md_nn" %in% names(obj@graphs)){
            obj <- Seurat::FindNeighbors(obj,
                                         reduction=reduction$md,
                                         dims=dims$md,
                                         k.param = k.param,
                                         graph.name ="md_nn")
        }
        if(!"ld_nn" %in% names(obj@graphs)){
            obj <- Seurat::FindNeighbors(obj,
                                         reduction=reduction$ld,
                                         dims=dims$ld,
                                         k.param = k.param,
                                         graph.name ="ld_nn")
        }
        md_nn <- obj@graphs$md_nn
        ld_nn <- obj@graphs$ld_nn 
        # run_knn_overlap_out$data$rnn$LD$graph$idx
        ### Iterate over each ID type 
        #RANN::nn2()
        consist_dt <- lapply(stats::setNames(id_types,id_types),
                             function(id_type){
                                 
         ids_dt <- data.table::as.data.table(
             obj@meta.data[,unique(c("orig.ident","id",id_type,"seurat_clusters"))]
         )[!is.na(get(id_type)),]
         dup_ids <- ids_dt[,c("n","n_clusters"):=list(.N,data.table::uniqueN(seurat_clusters)),
                           keyby=c(id_type)][n>1,]
         ids <- unique(ids_dt[get(id_type) %in% dup_ids[[id_type]]]$id)
         ## limit the number of IDs to avoid crashing Rstudio
         ids <- head(ids,max_ids)
         dup_ids <- dup_ids[id %in% ids]
         if(length(ids)==0){
             message(id_type," : ",length(ids)," duplicated IDs found. Reurning NULL.")
             return(NULL)
         } else{
             message(id_type," : ",length(ids)," duplicated IDs found.")
         }
         ## Subset obj 
         # obj_tmp <- obj[,unique(dup_ids$id)] 
         # obj_tmp <- Seurat::FindNeighbors(obj_tmp@assays[[1]],## causes immediate crashing: Dont run
         #                                    k.param = k.param,
         #                                    graph.name ="hd_nn")
         # obj_tmp <- Seurat::FindNeighbors(obj_tmp,
         #                                   reduction=reduction_md,
         #                                   dims=dims,
         #                                   k.param = k.param,
         #                                   graph.name ="md_nn")
         #  obj_tmp <- Seurat::FindNeighbors(obj_tmp,
         #                                   reduction=reduction_ld,
         #                                   dims=1:2,
         #                                   k.param = k.param,
         #                                   graph.name ="ld_nn")
         ## High-dimensional correlation
         # head(Seurat::VariableFeatures(obj),5000)
         dup_ids_cor_hd <- cor_fun(
             Seurat::GetAssayData(obj)[,ids,drop=FALSE])
         diag(dup_ids_cor_hd) <- NA
         ## Mid-dimensional correlation
         obsm_md <- obj@reductions[[reduction$md]]@cell.embeddings
         dup_ids_cor_md <- cor_fun(obsm_md[ids,,drop=FALSE], t=TRUE)
         diag(dup_ids_cor_md) <- NA
         ## Low-dimensional correlation
         obsm_ld <- obj@reductions[[reduction$ld]]@cell.embeddings
         dup_ids_cor_ld <- cor_fun(obsm_ld[ids,,drop=FALSE], t=TRUE)
         diag(dup_ids_cor_ld) <- NA
         ## Compute stats
         # md_agg <- orthogene:::aggregate_rows(X=as(md_nn,"sparseMatrix"),
         #                                      as_DelayedArray=FALSE,
         #                                      groupings = dup_ids$hp_id) 
         # md_agg <-sparseMatrixStats::rowAvgsPerColSet(as(md_nn,"sparseMatrix"),
         #                                              S = as.integer(dup_ids$hp_id))
         # ids <- intersect(colnames(run_knn_overlap_out$data$rnn$LD$data), 
         #                  dup_ids[hp_id==dup_ids$hp_id[1]]$id)
         # nn_hd <- rnndescent::rnnd_query(index = run_knn_overlap_out$data$rnn$LD$graph,
         #                                 query = run_knn_overlap_out$data$rnn$LD$data[,colinters]
         #                                 ) 
         message("Computing aggregate stats.")
         dup_ids_count <- dup_ids[,list(
             n_ids=unique(n),
             n_clusters=length(unique(seurat_clusters)),
             #nn_sim_hd= Matrix::mean(hd_nn[id[1],id[-1]]),
             nn_sim_md= Matrix::mean(md_nn[id[1],id[-1]]),
             nn_sim_ld= Matrix::mean(ld_nn[id[1],id[-1]]),
             cor_hd=Matrix::mean(dup_ids_cor_hd[id,id], na.rm = TRUE),
             cor_md=Matrix::mean(dup_ids_cor_md[id,id], na.rm = TRUE),
             cor_ld=Matrix::mean(dup_ids_cor_ld[id,id], na.rm = TRUE)
         ),
         by=c("id_type")]
         dup_ids_count[,c("baseline_cor_hd","baseline_cor_md","baseline_cor_ld",
                          "dim_hd","dim_md","dim_ld"):=list(
                              Matrix::mean(dup_ids_cor_hd, na.rm = TRUE),
                              Matrix::mean(dup_ids_cor_md, na.rm = TRUE),
                              Matrix::mean(dup_ids_cor_ld, na.rm = TRUE),
                              nrow(obj),
                              ncol(obsm_md),
                              ncol(obsm_ld)
                          )]
         dup_ids_count[,cluster_consistency_n1:=sum(n_clusters==1)/.N*100]  
         dup_ids_count[,cluster_consistency_n2:=sum(n_clusters<=2)/.N*100]  
         dup_ids_count[,cluster_consistency_n3:=sum(n_clusters<=3)/.N*100]  
         dup_ids_count|> data.table::setnames(id_type,"id_alt")
         gc()
         dup_ids_count
     }) |>data.table::rbindlist(idcol = "id_type")  
        ## Compute difference between baseline and identical ID correlation
        consist_dt[,diff_cor_hd:=cor_hd-baseline_cor_hd]
        consist_dt[,diff_cor_md:=cor_md-baseline_cor_md]
        consist_dt[,diff_cor_ld:=cor_ld-baseline_cor_ld]
        #### REWRITE SO ITS NOT HARD-CODED!
        nn_props <- 
            list(
                md=list(
                    all=round(sum(consist_dt$nn_sim_md>0)/
                                  nrow(consist_dt)*100,2),
                    HPO=round(sum(consist_dt[id_type=="hp_id"]$nn_sim_md>0)/
                                  nrow(consist_dt[id_type=="hp_id"])*100,2),
                    MONDO=round(sum(consist_dt[id_type=="mondo_id"]$nn_sim_md>0)/
                                    nrow(consist_dt[id_type=="mondo_id"])*100,2),
                    EFO=round(sum(consist_dt[id_type=="efo_id"]$nn_sim_md>0)/
                                  nrow(consist_dt[id_type=="efo_id"])*100,2)
                ),
                ld=list(
                    all=round(sum(consist_dt$nn_sim_ld>0)/
                                  nrow(consist_dt)*100,2),
                    HPO=round(sum(consist_dt[id_type=="hp_id"]$nn_sim_ld>0)/
                                  nrow(consist_dt[id_type=="hp_id"])*100,2),
                    MONDO=round(sum(consist_dt[id_type=="mondo_id"]$nn_sim_ld>0)/
                                    nrow(consist_dt[id_type=="mondo_id"])*100,2),
                    EFO=round(sum(consist_dt[id_type=="efo_id"]$nn_sim_ld>0)/
                                  nrow(consist_dt[id_type=="efo_id"])*100,2)
                ),
                random=round(Matrix::mean(Seurat::as.Graph(obj@graphs$md_nn), 
                                          na.rm = TRUE)*100, 2)
            )  
        #### Cache ####
        out <- list(consist_dt=consist_dt,
                    nn_props=nn_props)
        KGExplorer::cache_save(out, save_path = save_path)
        return(out)
    } 
}