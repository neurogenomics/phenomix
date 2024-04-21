#' Plot factors: Sankey
#' 
#' Create a Sankey plot from enrichment results performed on the 
#' feature weights and sample loadings from a factorization model.
#' @inheritParams ggsankey::geom_sankey
#' @inheritDotParams ggsankey::geom_sankey
#' @inheritParams ggplot2::scale_fill_discrete
#' @export
plot_factors_sankey <- function(factor.traits,
                                factor.celltypes,
                                varm = NULL, 
                                max_factors = NULL,
                                q_threshold=0.05,
                                logfc_threshold=NULL,
                                value_var="q",
                                drop=TRUE,
                                flow.alpha = .8,
                                label_size = 3,
                                label_color="white",
                                label_fill=ggplot2::alpha("black",.7),
                                xtext_size=12,
                                show.legend=FALSE,
                                ...
                                ){
    
    p_adjust_all <- factor_num <- NULL;
    
    X.traits <- data.table::dcast.data.table(factor.traits,
                                             formula = factor ~ name, 
                                             value.var = "p_adjust_all", 
                                             fun.aggregate = mean,
                                             fill=1)|> 
        KGExplorer::dt_to_matrix()
    hc.traits <- X.traits|>
        Matrix::t()|>
        stats::dist()|>
        stats::hclust()
    X.celltypes <- data.table::dcast.data.table(factor.celltypes,
                                                formula = factor ~  cl_name, 
                                                value.var = "q",
                                                fun.aggregate = mean,
                                                fill=1)|> 
        KGExplorer::dt_to_matrix()
    hc.celltypes <- X.celltypes|>
        Matrix::t() |>
        stats::dist()|>
        stats::hclust()
    if(!is.null(varm)){
        hc.factors <- varm |>
            Matrix::t() |>
            stats::dist()|>
            stats::hclust()
    } else {
        hc.factors <-
            Seurat::RowMergeSparseMatrices(X.traits,X.celltypes)|>
            stats::dist()|>
            stats::hclust()
    } 
    factor.annot <- merge(factor.traits[p_adjust_all<q_threshold,], 
                        factor.celltypes[q<q_threshold & mean_q<q_threshold,], 
                        by=c("factor","factor_num"), all.x = TRUE, 
                        allow.cartesian = TRUE) 
    MSTExplorer::add_logfc(factor.annot)
    if(!is.null(max_factors)){
        factor.annot <- factor.annot[factor_num<=max_factors]
    }
    if(!is.null(logfc_threshold)){
        factor.annot <- factor.annot[abs(logFC)>logfc_threshold &
                                     is.finite(logFC)]
    }
    dt <- ggsankey::make_long(factor.annot,
                              cl_name, factor_num, name, 
                              value = value_var)
    dt$value <- 1-dt$value 
    lvls <- c(stringr::str_split(hc.factors$labels[hc.factors$order],
                                 "_",simplify = TRUE)[,2],
              hc.traits$labels[hc.traits$order],
              hc.celltypes$labels[hc.celltypes$order]
    ) 
    dt$node <- factor(dt$node, levels = lvls, ordered = TRUE)
    dt$next_node <- factor(dt$next_node, levels = lvls, ordered = TRUE)
    #### Create plot ####
    gg_sankey <- ggplot2::ggplot(dt,
                    ggplot2::aes(x = x,
                                 next_x = next_x,
                                 node = node,
                                 next_node = next_node,
                                 fill = factor(node),
                                 label = node)
    ) +
        ggsankey::geom_sankey(show.legend = show.legend,
                              flow.alpha = flow.alpha,
                              ...) +
        ggsankey::geom_sankey_label(size = label_size, 
                                    color = label_color, 
                                    fill = label_fill) +
        ggplot2::scale_fill_discrete(drop=drop) +
        ggplot2::labs(x=NULL) +
        ggplot2::scale_x_discrete(labels=c(
            expression("Cell types   " %<-% "   factor gene weights   " %<-% ""),
            "Factor",
            expression("" %->% "   factor trait loadings   " %->% "   Traits" ))
        ) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(hjust = c(0,.5,1),
                                                           size = xtext_size)) +
        theme_nightlight()
    #### Return ####
    return(list(
        data=dt,
        plot=gg_sankey
    ))
    
    # meta <- merge(
    #   data.table::data.table(obj@meta.data|>dplyr::select(!id),
    #                          keep.rownames = "id"),
    #   scKirby::get_obsm(obj, keys = reduction_md, n=1)^2|>
    #     data.table::as.data.table(keep.rownames = "id")|>
    #     data.table::melt.data.table(
    #       id.vars = "id",
    #       variable.name = "factor",
    #       value.name = "loading"),
    #   by="id"
    # )
    
    # id_to_cluster <- tidygraph::as_tbl_graph(
    #     data.table::data.table(colnames(obj),
    #                            cluster=paste0("cluster",obj$seurat_clusters),
    #                            value=1,
    #                            edge_type="id_to_cluster")
    # )
    # id_to_factor <- tidygraph::as_tbl_graph(
    #     (scKirby::get_obsm(obj, keys = reduction_md, n=1)^2|>
    #          data.table::as.data.table(keep.rownames = "id")|>
    #          data.table::melt.data.table(
    #              id.vars = "id",
    #              variable.name = "factor",
    #              value.name = "value")
    #     )[,edge_type:="id_to_factor"]
    # )
    # factor_to_celltype <- tidygraph::as_tbl_graph(
    #     factor.celltypes[q<.05,c("factor","cl_name","p")][,edge_type:="factor_to_celltype"] |>
    #         data.table::setnames("p","value")
    # )
    # 
    # visNetwork::visIgraph(factor_to_celltype) |>
    #     # visNetwork::visHierarchicalLayout()|> 
    #     visNetwork::visIgraphLayout(layout = "layout_with_sugiyama")|>
    #     visNetwork::visNodes() |>
    #     visNetwork::visEdges(color = list(opacity=.1))|>
    #     visNetwork::visInteraction(hover = TRUE)
    # 
    # g <- tidygraph::graph_join(id_to_cluster,
    #                            id_to_factor)|>
    #     tidygraph::graph_join(factor_to_celltype)
    # # g <- tidygraph::sample_n(g, 1000)
    # g.dt <- KGExplorer::graph_to_dt(g)
    # nodes <- KGExplorer::graph_to_dt(g, what = "nodes")
    # edges <- KGExplorer::graph_to_dt(g, what = "edges")
    # edges[,c("from","to"):=list(from-1,to-1)] # networkD3 requirees 0-indexing
    # p <- networkD3::sankeyNetwork(Nodes = nodes,
    #                               Links = edges, 
    #                               Source = "from",
    #                               Target = "to", 
    #                               Value = "weight", 
    #                               NodeID = "name")
    # p
    # 
    # 
    # plotly::plot_ly(type="sankey",g)
    
    # X <- data.table::dcast.data.table(data = factor.celltypes, 
    #                                   formula = factor~cl_name,
    #                                   value.var = "p",
    #                                   fun.aggregate = mean, na.rm = TRUE)|>
    #   KGExplorer::dt_to_matrix()
    # 
    # heatmaply::heatmaply(as.matrix(X))
    
}