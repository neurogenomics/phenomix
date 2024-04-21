run_ontological_velocity_sampled <- function(obj,
                                             ont,
                                             iterations,
                                             graph=tail(obj@graphs,1)[[1]],
                                             k=20,
                                             min_dist=0.01,
                                             max_dist,
                                             seed = sample(colnames(obj),1),
                                             nn = 
                                                 tail(
                                                     sort(
                                                         graph[,seed][graph[,seed]>min_dist]
                                                     ),
                                                     k),
                                             agg_fun=stats::weighted.mean
                ){  
    run_ontological_velocity_sampled_i <- function(...){
        id1 <- id2 <- NULL;
        # messager("seed=",shQuote(seed))
        nns <- names(nn) 
        if(length(nns)<2) return(NULL)
        emb <- (Seurat::Embeddings(obj, reduction="umap")|>data.frame())[nns,]
        #### remove suffixes from make.unique ####
        nns <- map_id_sep(
            stringr::str_split(nns,"[.]", simplify = TRUE)[,1]
        ) 
        nns <- intersect(nns,ont@terms)
        if(length(nns)<2) return(NULL)
        emb$id <- map_id_sep(
            stringr::str_split(rownames(emb),"[.]", simplify = TRUE)[,1]
        ) 
        combos <- expand.grid(id1=nns,
                              id2=nns) |>
            data.table::data.table()
        combos <- combos[id1!=id2,]
        if(nrow(combos)==0) return(NULL)
        d <- suppressMessages(
            simona::longest_distances_via_LCA(ont,
                                              terms = nns)    
        )  
        d2 <- suppressMessages(
            simona::shortest_distances_directed(ont,
                                                terms = nns)
        )
        get_dist <- function(id1,id2,ont,d){
            id1 <- as.character(id1)
            id2 <- as.character(id2)
            d[id1,id2]
            ## *(if(d2[id1,id2]<0) -1 else 1)
        }
        combos[,ont_distance:=get_dist(id1,id2,ont,d), by=.I]
        combos[,ont_direction:=ifelse(get_dist(id1,id2,ont,d2)>0,1,-1), by=.I]
        tmp<- merge(combos,
                    emb,
                    by.x="id1",
                    by.y="id")|>
            merge(emb,
                  by.x="id2",
                  by.y="id",
                  suffixes=c(".id1",".id2"))
        tmp[,knn_sim:=graph[id1,id2], by=.I]
        ## (y₂ - y₁)/(x₂ - x₁)
        tmp[,dx:=ifelse(ont_direction>0,
                        umap_1.id1-umap_1.id2,
                        umap_1.id2-umap_1.id1)]
        tmp[,dy:=ifelse(ont_direction>0,
                        umap_2.id1-umap_2.id2,
                        umap_2.id2-umap_2.id1)]
        tmp[,slope:=dy/dx]
        tmp[,ont_distance_scaled:=ont_distance/max_dist]
        tmp[,weight:=ont_distance_scaled/(1-knn_sim)]
        tmp[,n_points:=length(nns)]
        # tmp
        ### aggregate
        if(!is.null(agg_fun)){
            cols <- sapply(tmp, is.numeric)
            cols <- names(cols)[cols]
            tmp_mean <- tmp[, lapply(.SD, agg_fun,weight), .SDcols = cols]
            return(tmp_mean)   
        }else{
            return(tmp)
        }
    }
    BPPARAM <- KGExplorer::set_cores() 
    slopes <- BiocParallel::bplapply(BPPARAM = BPPARAM,
                                     seq(iterations), 
                                     function(i){
         messager("iteration",i,v=!BPPARAM$progressbar)
         tmp <- run_ontological_velocity_sampled_i(obj = obj,
                                                   ont = ont, 
                                                   k = k,
                                                   graph = graph,
                                                   agg_fun = agg_fun,
                                                   max_dist = max_dist)
         # if(!is.null(tmp)) tmp[,i:=i]
         tmp
         }) |> data.table::rbindlist()
    return(slopes)
} 
