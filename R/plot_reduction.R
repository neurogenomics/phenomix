#' Plot reduction
#'
#' Plot dimensionality reduction results
#' from \code{run_<method>} or a Seurat object.
#'
#' @param metadata Metadata associated with \code{obj}. 
#' Will be automatically extracted from \code{obj} if it is a
#'  \code{Seurat} object.
#' @param fix_rownames Replace non-ASCII characters in row names.
#' @param color_var \code{metadata} variable to use as point color.
#' @param label_var \code{metadata} variable to use as point labels
#' @param size_var \code{metadata} variable to use as point size.
#' @param labels Whether to show labels using \link[ggrepel]{geom_text_repel}.
#' @param show_plot Whether to print the plot.
#' @param point_alpha Opacity of data points.
#' 
#' @inheritParams extract_embeddings
#' @inheritParams ggrepel::geom_text_repel
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
                           x_dim = 1,
                           y_dim = 2,
                           fix_rownames = FALSE,
                           color_var = NULL,
                           label_var = NULL,
                           size_var = NULL,
                           labels = !is.null(label_var),
                           max.overlaps = 30,
                           show_plot = TRUE,
                           point_alpha = .6) { 
    if (is.null(metadata)) {
        metadata <- extract_metadata(obj = obj)
    }
    if (is.null(rownames(metadata))) {
        stop("metadata must have rownames to merge with obj.")
    } 
    #### Prepare data ####
    embeddings <- extract_embeddings(
        obj = obj,
        reduction = reduction
    ) 
    if(x_dim>ncol(embeddings) || x_dim<1) {
        stop("x_dim must be an integer >= the number of embeddings dimensions.") 
    }
    if(y_dim>ncol(embeddings) || y_dim<1) {
        stop("y_dim must be an integer >= the number of embeddings dimensions.") 
    }
    #### Fix rownames to match metadata rownames ####
    if (fix_rownames) {
        metadata <- metadata %>%
            `rownames<-`(gsub("-|[:]|[.]", "_", rownames(metadata)))
        embeddings <- embeddings %>%
            `rownames<-`(gsub("-|[:]|[.]", "_", rownames(embeddings)))
    }
    plot_df <- merge(
        x = embeddings,
        y = metadata[,!startsWith(colnames(metadata),"UMAP")],
        by = 0
    )
    if(nrow(plot_df)==0){
        stop("No rows in merged metadata + embeddings data.")
    } else {
        messager("Constructed plot_df:",
                 formatC(nrow(plot_df),big.mark = ","),"rows x",
                 formatC(ncol(plot_df),big.mark = ","),"columns"
                 )
    }

    gp <- ggplot(plot_df, aes_string(
        x = colnames(embeddings)[x_dim],
        y = colnames(embeddings)[y_dim],
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
