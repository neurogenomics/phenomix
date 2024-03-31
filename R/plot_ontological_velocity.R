#' Plot ontological velocity
#' 
#' Plot ontological velocity.
#' 
#' @import data.table
#' @export
#' @examples
#' obj <- get_HPO()
#' ont <- KGExplorer::get_ontology("hp")
#' out <- plot_ontological_velocity(obj,ont)
plot_ontological_velocity <- function(obj,
                                      ont=KGExplorer::get_ontology("hp"),
                                      id_col=NULL,
                                      iterations=5000,
                                      k=30,
                                      agg_fun=stats::weighted.mean,
                                      show_plot=TRUE){
    # https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html
    # https://ouyanglab.com/singlecell/clust.html
    # https://www.bioconductor.org/packages/devel/bioc/vignettes/Nebulosa/inst/doc/nebulosa_seurat.html
    # https://eliocamp.github.io/metR/reference/geom_arrow.html
    # brew install boost
    # obj2 <- SeuratWrappers::RunVelocity(obj)
    requireNamespace("ggplot2")
    requireNamespace("metR")
    
    g <- KGExplorer::ontology_to(ont,"tbl_graph") 
    ## Step 1. Sample K-nearest neighbours for 100 samples 
    ## of the low-dimensional graph embedding (UMAP).
    if(!is.null(id_col)){
        if(id_col=="hp_id"){
            obj$hp_id <- data.table::fcoalesce(
                obj$hp_id,
                ifelse(grepl("^HP",colnames(obj)),
                       colnames(obj),NA))   
        }
        obj <- obj[,!is.na(obj@meta.data[[id_col]])]
        colnames(obj) <- make.unique(replace_char_fun(obj$hp_id))
    }  
    dist_lca <- simona::longest_distances_via_LCA(ont,
                                                  terms = ont@terms)
    # dist_nca <- simona::shortest_distances_via_NCA(ont,
    #                                                terms = ont@terms)
    max_dist <- max(dist_lca) 
    compute_velocity <- function(obj,
                                 ont,
                                 graph=tail(obj@graphs,1)[[1]],
                                 k=30,
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
        id1 <- id2 <- NULL;
        # messager("seed=",shQuote(seed))
        nns <- names(nn) 
        if(length(nns)<2) return(NULL)
        emb <- (obj@reductions$umap@cell.embeddings |>data.frame())[nns,]
        #### remove suffixes from make.unique ####
        nns <- replace_char_fun(
            stringr::str_split(nns,"[.]", simplify = TRUE)[,1]
        ) 
        nns <- intersect(nns,ont@terms)
        if(length(nns)<2) return(NULL)
        emb$id <- replace_char_fun(
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
    graph <- tail(obj@graphs,1)[[1]]
    BPPARAM <- KGExplorer::set_cores() 
    slopes <- BiocParallel::bplapply(BPPARAM = BPPARAM,
                                     seq(iterations), 
                                     function(i){
        messager("iteration",i,v=!BPPARAM$progressbar)
        tmp <- compute_velocity(obj = obj,
                                ont = ont, 
                                k = k,
                                graph = graph,
                                agg_fun = agg_fun,
                                max_dist = max_dist)
        # if(!is.null(tmp)) tmp[,i:=i]
        tmp
    }) |> data.table::rbindlist()
    
    slopes <- slopes[!is.na(slope),]
    ### aggregate
    cols <- sapply(slopes, is.numeric)
    cols <- names(cols)[cols]
    # slopes_mean <- slopes[, lapply(.SD, stats::weighted.mean,weight, na.rm=TRUE), 
    #                       .SDcols = cols]
    slopes[,slope_capped:=pmin(pmax(log10(abs(slope)),-3),3)]
    slopes[,slope_logged:=log10(abs(slope))]
    
    
    # slopes_melt <- data.table::melt(slopes,
    #                                 id.vars = c("id1","id2","slope","weight"),
    #                                 measure.vars = list(
    #                                     umap_1=c("umap_1.id1","umap_1.id2"),
    #                                     umap_2=c("umap_2.id1","umap_2.id2")
    #                                                     )
    #                                 )
   
    p <- ggplot2::ggplot() +
        # ggplot2::geom_point(data=slopes_melt[!is.na(weight) & weight>4],#obj@reductions$umap@cell.embeddings, 
        #                     ggplot2::aes(x = umap_1,
        #                                  y = umap_2,
        #                                  color=weight),
        #                     alpha=.5) +
        
        metR::geom_arrow(data=slopes,
                         pivot = 1,# 0 beginning, 1=end
                         # min.mag = 0,
                         # arrow.length=1,
                         alpha=.7,
                         ggplot2::aes(x=umap_1.id1, y=umap_2.id1,
                                      dx = dx, dy = dy,
                                      alpha=1-ont_distance_scaled,
                                      color=slope
                                      )) +
        # metR::scale_mag(max_size = 1) +
        # paletteer::scale_color_paletteer_c(palette = "pals::kovesi.diverging_bwr_40_95_c42") +
        ggplot2::scale_color_viridis_c(option = "plasma",alpha = .5)+ 
        ggplot2::theme_bw()
    if(show_plot) methods::show(p)
    return(list(
        plot=p,
        data=slopes
    ))
    
    # Seurat::FeaturePlot(obj, pt.size = 1,
    #                     # cols=c("grey","blue","red"),
    #                     features  = "ontLvl") +
    #     ggplot2::scale_color_viridis_c(option="plasma")
    
    # ggplot2::ggplot() +
    # ggplot2::geom_point(data=obj@reductions$umap@cell.embeddings,
    #                     ggplot2::aes(
    #     x = umap_1,
    #     y = umap_2),
    #     alpha=.1) +
    # ggplot2::geom_curve(data=slopes[abs(dx)>0.001 & abs(dy)>0.001,],
    #                       ggplot2::aes(
    #     x = umap_1.id1,
    #     y = umap_2.id1,
    #     xend = umap_1.id2,
    #     yend = umap_2.id2,
    #     color = abs(shortest_distance),
    #     alpha=abs(shortest_distance),
    #     linewidth = abs(shortest_distance)
    #     ),
    #     arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "inches"))
    #     )
        
    
    # lvls <- seq(max(ont@elementMetadata$ontLvl))
    # for (l in lvls){
    #     ont <- KGExplorer::add_ancestors(ont, 
    #                                      prefix = paste0("ancestor",l),
    #                                      lvl=l, 
    #                                      force_new=TRUE)
    # }
    # 
    # meta <- merge(
    #       obj@meta.data,
    #       ont@elementMetadata[,],
    #       by.x=0,
    #       by.y="id",
    #       all.x = TRUE
    # )
    # rownames(meta) <- meta$Row.names
    # obj@meta.data <- meta 
    # clustree::clustree(obj,
    #                    prefix = "ancestor") 
    
}
