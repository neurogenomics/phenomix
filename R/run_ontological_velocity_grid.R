run_ontological_velocity_grid <- function(obj,
                                          ont, 
                                          graph_name,
                                          max_dist=0.2,
                                          iterations=100,
                                          agg_fun=function(x){mean(x,na.rm=TRUE)},
                                          reduction="umap"
                                          ){  
    requireNamespace("proxy")
    id1 <- id2 <- NULL;
    graph <- obj@graphs[[graph_name]]
    dist_lca <- simona::longest_distances_via_LCA(ont,
                                                  terms = ont@terms)
    # dist_nca <- simona::shortest_distances_via_NCA(ont,
    #                                                terms = ont@terms)
    dist_lca.max <- max(dist_lca)  
    embedding <- Seurat::Embeddings(obj,
                                    reduction=reduction)|>data.frame() 
    colnames(embedding) <- c("x","y")  
    fix_id <- function(x){
        as.character(map_id_sep(gsub("*\\.","",x)))
    }
    embedding$id <- rownames(embedding)
    embedding$id_mapped <- fix_id(rownames(embedding))
    grid_dt <- expand.grid(
        x=seq(min(embedding[[1]]), max(embedding[[1]]),
              length.out = round(sqrt(iterations))),
        y=seq(min(embedding[[2]]), max(embedding[[2]]),
              length.out = round(sqrt(iterations)))
    ) |>data.table::data.table()
    ## Compute euclidean distance between xy in df and all points in embedding
    get_dist <- function(x,y,embedding){
        sqrt((x-embedding[[1]])^2 + (y-embedding[[2]])^2)
    }
    grid_matches <- grid_dt[,list(
        x=x,
        y=y,
        id=rownames(embedding)[get_dist(x,y,embedding)<=max_dist],
        dist=get_dist(x,y,embedding)[get_dist(x,y,embedding)<=max_dist]
        ), 
        by=.I]
    grid_matches[,id_mapped:=fix_id(id)]
    grid_matches[id_mapped %in% ont@terms, n:=data.table::uniqueN(id),by=c("I")] 
    ## Only run for grid points with more than one data point 
    ## (so distances can be computed)
    # nns <- intersect(nns,ont@terms)
    grid_ids <- unique(grid_matches[n>1]$I)
    messager("Running ontological velocity analysis with",
             length(grid_ids),"iterations.")
    BPPARAM <- KGExplorer::set_cores() 
    slopes <- BiocParallel::bplapply(BPPARAM = BPPARAM,
                                     grid_ids,
                                     function(i){              
        gm <- grid_matches[I==i,][id_mapped %in% ont@terms,]
        id_mapped <- unique(gm$id_mapped)
        ## Create all pairwise combinations between nns
        combos <- expand.grid(id1=gm$id,
                              id2=gm$id) |>
            data.table::data.table()
        combos <- combos[id1!=id2,]
        sim <- suppressMessages(
            simona::term_sim(dag = ont,
                             terms = id_mapped)  
        )   
        dist_lca <- suppressMessages(
            simona::longest_distances_via_LCA(ont,
                                              terms = id_mapped)    
        )  
        dist_direction <- suppressMessages(
            simona::shortest_distances_directed(ont,
                                                terms = id_mapped)
        )
        get_metric <- function(id1,id2,x){
            x[fix_id(id1),fix_id(id2)]
        }
        combos[,ont_sim:=get_metric(id1,id2,sim), by=.I]
        combos[,ont_dist:=get_metric(id1,id2,dist_lca), by=.I]
        combos[,ont_direction:=ifelse(get_metric(id1,id2,dist_direction)>0,1,-1), by=.I]
        tmp<- merge(combos,
                    embedding,
                    by.x="id1",
                    by.y="id")|>
            merge(embedding,
                  by.x="id2",
                  by.y="id",
                  suffixes=c(".id1",".id2"))
        tmp[,n_ids:=length(id_mapped)]
        embed_dist_X <- proxy::dist(embedding[embedding$id %in% gm$id,c(1,2)])|>as.matrix()
        embed_sim_X <- proxy::simil(embedding[embedding$id %in% gm$id,c(1,2)])|>as.matrix()
        tmp[,graph_dist:=graph[id1,id2], by=.I]
        tmp[,embed_dist:=embed_dist_X[id1,id2], by=.I]
        tmp[,embed_sim:=embed_sim_X[id1,id2], by=.I]
        ## (y₂ - y₁)/(x₂ - x₁)
        tmp[,dx:=ifelse(ont_direction>0,
                        x.id1-x.id2,
                        x.id2-x.id1)]
        tmp[,dy:=ifelse(ont_direction>0,
                        y.id1-y.id2,
                        y.id2-y.id1)]
        ## Compute euclidean distance
        tmp[,slope:=dy/dx][is.infinite(slope),slope:=NA] 
        tmp[,ont_dist_scaled:=ont_dist/dist_lca.max]
        tmp[,ont_velocity:=(embed_sim/max(ont_sim,.Machine$double.xmin))-1]
        if(nrow(tmp[is.infinite(ont_velocity)])>0){
            messager(nrow(tmp[is.infinite(ont_velocity)]),
                     "infinite values in ont_velocity, setting to NA.")
            tmp[is.infinite(ont_velocity),ont_velocity:=NA]
        }
        ### Add extra info 
        tmp <- cbind(I=i,
                     max_dist.grid=max_dist,
                     mean_dist.grid=mean(gm$dist),
                     x=gm$x[1],
                     y=gm$y[1],
                     tmp)
        ### aggregate
        if(!is.null(agg_fun)){
            cols <- sapply(tmp, is.numeric)
            cols <- names(cols)[cols]
            if(identical(agg_fun,stats::weighted.mean)){
                tmp_mean <- tmp[, lapply(.SD, agg_fun,ont_velocity), .SDcols = cols]    
            }else {
                tmp_mean <- tmp[, lapply(.SD, agg_fun), .SDcols = cols]
            } 
            return(tmp_mean)   
        }else{
            return(tmp)
        }
    }) |> data.table::rbindlist(fill=TRUE)
    #### Add grid points with no surrounding data points 
    slopes <- rbind(
        slopes,
        grid_matches[!I %in% slopes$I,],
        fill=TRUE
    )
    return(slopes)
} 
