#' Plot reduction
#'
#' Plot dimensionality reduction results
#' from \code{run_<method>} or a Seurat object.
#'
#' @inheritParams extract_embeddings
#'
#' @return ggplot object
#'
#' @export
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr %>%
plot_reduction <- function(obj,
                           reduction = NULL,
                           metadata = NULL,
                           fix_rownames = FALSE,
                           color_var = NULL,
                           label_var = NULL,
                           size_var = NULL,
                           labels = TRUE,
                           max.overlaps = 30,
                           show_plot = TRUE,
                           point_alpha = .6) {
    if (is.null(label_var)) labels <- FALSE
    if (is.null(metadata)) {
        metadata <- extract_metadata(obj = obj)
    }
    if (is.null(rownames(metadata))) {
        stop("metadata must have rownames to merge with obj.")
    }
    #### Fix rownames to match metadata rownames ####
    if (fix_rownames) {
        metadata <- metadata %>%
            `rownames<-`(gsub("-|[:]|[.]", "_", rownames(metadata)))
    }
    #### Prepare data ####
    embeddings <- extract_embeddings(
        obj = obj,
        reduction = reduction
    )
    umap_df <- merge(
        x = embeddings,
        y = metadata,
        by = 0
    )

    gp <- ggplot(umap_df, aes_string(
        x = colnames(embeddings)[1],
        y = colnames(embeddings)[2],
        color = color_var, label = label_var
    )) +
        geom_point(aes_string(size = size_var), alpha = point_alpha) +
        theme_bw()
    if (labels) {
        gp <- gp + ggrepel::geom_text_repel(max.overlaps = max.overlaps)
    }
    if (show_plot) print(gp)
    return(gp)
}
