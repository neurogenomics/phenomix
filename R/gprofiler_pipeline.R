#' gprofiler_pipeline
#' 
#' @keywords internal
#' @importFrom dplyr arrange slice_head
gprofiler_pipeline <- function(markers_lab,
                               p_val_adj_thresh = .05,
                               markers_per_cluster = 50) {
    cluster <- p_val_adj <- avg_log2FC <- NULL;
    gp_list <- lapply(unique(markers_lab$cluster), function(x) {
        markers_sub <- subset(markers_lab, 
                              cluster == x & 
                                  p_val_adj < p_val_adj_thresh) |>
            dplyr::arrange(desc(avg_log2FC)) |>
            dplyr::slice_head(n = markers_per_cluster)
        g <- list(markers_sub$gene)
        names(g) <- paste(unique(markers_sub$cluster),
            unique(markers_sub$enriched_words),
            sep = "_"
        )
        messager("cluster = ", x, "; genes = ",
                length(markers_sub$gene))
        gp_res <- tryCatch(expr = {
            gprofiler2::gost(query = g)
        }, error = function(e) NULL)
        return(gp_res)
    }) |> `names<-`(unique(markers_lab$enriched_words))
    return(gp_list)
}
