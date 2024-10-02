#' Run K-nearest neighbors overlap
#' 
#' Computes KNN graphs based on different data types:
#' \itemize{
#' \item{ONT: }{Ontological similarity}
#' \item{HD: }{High-dimensional data (e.g. gene associations)}
#' \item{MD: }{Mid-dimensional data (e.g. PCA factors)}
#' \item{LD: }{Low-dimensional data (e.g. UMAP components)}
#' \item{Random: }{Random data to compare to as a baseline}
#' }
#' It then computes the overlap between the KNN graphs and returns the
#'  distributions as a data table.
#' Finally, it plots the overlap results.
#' @param save_full Save the full results. 
#' NOTE: this can be very large (multiple Gb) and take a long time to save.
#' @inheritParams rnndescent::rnnd_build
#' @inheritDotParams rnndescent::rnnd_build
#' @export
run_knn_overlap <- function(obj,
                            ont,
                            id_col,
                            assay=Seurat::DefaultAssay(obj),
                            k = 30,
                            nfeatures=NULL,
                            reduction_md="pca",
                            reduction_ld="umap",
                            reference="ONT",
                            n_threads = NULL,
                            low_memory=FALSE,  
                            show_plot = TRUE,
                            title=paste0("KNN overlap with ",reference),
                            verbose = TRUE,
                            save_full=FALSE,
                            save_path=NULL,
                            force_new=FALSE,
                            ...){
    requireNamespace("ggplot2")
    requireNamespace("rnndescent")
    id1 <- id2 <- overlap <- logFC <- NULL;
     
    #### Load existing results ####
    if(!is.null(save_path) &&
       file.exists(save_path) && 
       isFALSE(force_new)){
        message("Loading existing results from ",save_path)
        data <- readRDS(save_path) 
    } else {
    #### Run new analyses ####
        if(is.null(n_threads)) n_threads <- KGExplorer::set_cores()$workers
        
        rnn <- list()
        ## Ont
        top_ic <- rownames(ont@elementMetadata)[order(ont@elementMetadata$IC, 
                                                      decreasing=TRUE)]
        ontsim <- seurat_to_ontological_similarity(obj=obj,
                                                   id_col = id_col,
                                                   ont = ont,
                                                   ancestors = head(top_ic,Inf), 
                                                   return_assay=FALSE
        ) 
        rnn[["ONT"]] <- rnndescent::rnnd_build(Matrix::t(ontsim$Xsim),
                                               k = k,
                                               low_memory=low_memory,  
                                               n_threads=n_threads,
                                               verbose=verbose,
                                               ...)
        ids <- colnames(rnn[["ONT"]]$data)
        ## HD
        if(!is.null(assay)){
            hd <- Seurat::GetAssayData(obj, assay=assay)[,ids] 
            if(!is.null(nfeatures)){
                hd <- hd[Seurat::VariableFeatures(obj)[seq(nfeatures)],]
            }
            rnn[["HD"]] <- rnndescent::rnnd_build(Matrix::t(hd), 
                                                  k = k, 
                                                  n_threads=n_threads,
                                                  low_memory=low_memory, 
                                                  verbose = verbose,
                                                  ...)
            remove(hd); gc()
        } 
        ## MD
        if(!is.null(reduction_md)){ 
            md <- scKirby::get_obsm(obj, keys = reduction_md, n = 1)[ids,]
            rnn[["MD"]] <- rnndescent::rnnd_build(md,
                                                  k = k, 
                                                  n_threads=n_threads,
                                                  low_memory=low_memory,
                                                  verbose = verbose,
                                                  ...)
        }
        ## LD
        if(!is.null(reduction_ld)){ 
            ld <- scKirby::get_obsm(obj, keys = reduction_ld, n = 1)[ids,]
            rnn[["LD"]]  <- rnndescent::rnnd_build(ld, 
                                                   k = k,  
                                                   n_threads=n_threads,
                                                   low_memory=low_memory,
                                                   verbose = verbose,
                                                   ...)
        } 
        ## generate random neighbors
        M <- matrix(data = rnorm(length(ids)*2), 
                    nrow = length(ids),
                    ncol = 2,
                    dimnames = list(ids,NULL))
        rnn[["Random"]] <- rnndescent::rnnd_build(M,
                                                  k = k,
                                                  n_threads=n_threads,
                                                  low_memory=low_memory,
                                                  verbose = verbose,
                                                  ...)
        ### compute overlap for each pairwise combinations of graphs
        messager("Computing iterative neighbour overlap.")
        cmb <- combn(names(rnn),2)
        rnn_res <- lapply(seq(ncol(cmb)), function(i){
            o <- rnndescent::neighbor_overlap(idx1 = rnn[[cmb[1,i]]]$graph,
                                              idx2 = rnn[[cmb[2,i]]]$graph, 
                                              ret_vec = TRUE)
            data.table::data.table(
                id1=cmb[1,i],
                id2=cmb[2,i],
                overlap=o$overlap
            )
        })|>data.table::rbindlist()
        ## postprocess
        rnn_res[,id2:=factor(id2,levels=c("ONT","Random","HD","MD","LD"), 
                             ordered=TRUE)]
        rnn_res[id1=="ONT",logFC:=log2(overlap/rnn_res[id2=="Random",overlap])]
        rnn_res[is.infinite(logFC),logFC:=0]
        data <- list(rnn_res=rnn_res,
                     ids=ids)
        if(save_full | is.null(save_path)){
            data[["ontsim"]] <- ontsim
            data[["rnn"]] <- rnn
        }
        ### Save only the essentials 
        KGExplorer::cache_save(data, save_path)
        ### Add full data after saving is done
        data[["ontsim"]] <- ontsim
        data[["rnn"]] <- rnn
        remove(ontsim, rnn, rnn_res)
    } 
    ### Agg data
    data[["rnn_mean"]] <- data$rnn_res[id1==reference,
                                       lapply(.SD,mean, na.rm=TRUE),
                                       by=c("id2"),
                                       .SDcols = c("overlap","logFC")]
    rnn_res_agg <- data$rnn_res[,lapply(.SD,mean,na.rm=TRUE),
                                by=c("id1","id2"), 
                                .SDcols=c("overlap","logFC")]
    rnn_res_agg <- rbind(rnn_res_agg,
          rnn_res_agg[,.(id1=id2,id2=id1,overlap,logFC)])
    #### Plot ####
    ## bar plots
    gg_bar <- ggplot2::ggplot(data$rnn_res[id1==reference,], 
                              ggplot2::aes(x=id2, y=logFC, 
                                           fill=id2)) +
        ggplot2::labs(title=title,
                      subtitle=paste0("k=",k),
                      y="log2(fold-change) of KNN overlap",
                      x="Data Representation",
                      fill=NULL) +
        ggplot2::geom_boxplot() + 
        ggplot2::geom_violin() +
        ggplot2::geom_label(data = data$rnn_mean,
                            hjust=0,
                            fill= ggplot2::alpha("white",.75),
                            ggplot2::aes(
                                y=max(data$rnn_res[id1==reference,]$logFC)*.95,
                                label=paste0(
                                    "Averages:\n",
                                    "log(FC)=",round(logFC,1),"\n",
                                    "KNN overlap=",round(overlap*100,1),"%"
                                )
                            )
        ) +
        ggplot2::scale_fill_viridis_d(begin = .3) +
        ggplot2::theme_bw()
    
    ## histograms 
    gg_hist <- ggplot2::ggplot(data$rnn_res[id1==reference,], 
                               ggplot2::aes(x=overlap,
                                            fill=id2)) +
        ggplot2::facet_wrap(~id2,ncol = 1) +
        ggplot2::geom_histogram() +
        ggplot2::labs(title=title,
                      subtitle=paste0("k=",k),
                      y=paste0("KNN overlap with ",reference)) +
        ggplot2::scale_fill_viridis_d(begin = .3) +
        ggplot2::theme_bw() 
    ## tiles 
    gg_tile <- ggplot2::ggplot(rnn_res_agg) +
        ggplot2::geom_tile(ggplot2::aes(x=id1, y=id2, fill=overlap)) +
        ggplot2::labs(fill="KNN\noverlap") +
        ggplot2::theme_bw() +
        ggplot2::labs(x="Data Representation 1",
                      y="Data Representation 2") +
        ggplot2::scale_fill_viridis_c(option = "plasma") 
    
    return(list(
        data=data,
        plot=list(
            bar=gg_bar,
            hist=gg_hist,
            tile=gg_tile
        )
    ))
}


### Compute similarity and get infer KNN (inefficient but it ensures you have comparable KNN graphs)
# gcompare <- igraph::compare(g_md, g_ld)

### Literally ontological distance (n steps in graph)
# id_col="mondo_id"
# ont <- mondo
# ontDist <- simona::longest_distances_via_LCA(mondo, 
#                                               terms = intersect(mondo@terms,
#                                                                 obj$mondo_id))
# matched_meta <- obj@meta.data[
#         simona::dag_has_terms(ont,obj@meta.data[[id_col]]),
#     ]
# og_ids <- stats::setNames(matched_meta[[id_col]],
#                           rownames(matched_meta)
#   )[matched_meta[[id_col]] %in% colnames(ontDist)]
# ontDist <- ontDist[unname(og_ids),unname(og_ids)]
# og_ids_rev <- KGExplorer:::invert_dict(og_ids)
# colnames(ontDist) <- og_ids_rev[colnames(ontDist)]
# rownames(ontDist) <- og_ids_rev[rownames(ontDist)] 
# ont_rnn <- list(idx=, graph=)

