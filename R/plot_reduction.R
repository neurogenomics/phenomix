#' Plot reduction
#'
#' Plot dimensionality reduction results
#' from \code{run_<method>} or a Seurat object.
#'
#' @param obs Phenotype metadata associated with \code{obj}. 
#' Will be automatically extracted from \code{obj} if it is a
#'  \code{Seurat} object.
#' @param fix_rownames Replace non-ASCII characters in row names.
#' @param color_var \code{obs} variable to use as point color.
#' @param label_var \code{obs} variable to use as point labels
#' @param size_var \code{obs} variable to use as point size.
#' @param labels Whether to show labels using \link[ggrepel]{geom_text_repel}.
#' @param show_plot Whether to print the plot.
#' @param x_dim Which dimension to plot on the x-axis.
#' @param y_dim Which dimension to plot on the y-axis.
#' @param point_alpha Opacity of data points.
#' @inheritParams scKirby::get_obsm
#' @inheritParams ggrepel::geom_text_repel
#' @returns ggplot object
#'
#' @export 
#' @importFrom ggrepel geom_text_repel
#' @examples
#' obj <- get_HPO()
#' gp <- plot_reduction(obj = obj)
plot_reduction <- function(obj,
                           keys = NULL,
                           obs = NULL,
                           x_dim = 1,
                           y_dim = 2,
                           fix_rownames = FALSE,
                           color_var = NULL,
                           label_var = NULL,
                           size_var = NULL,
                           labels = !is.null(label_var),
                           max.overlaps = 30,
                           show_plot = TRUE,
                           point_alpha = .6,
                           verbose = TRUE) { 
    # devoptera::args2vars(plot_reduction)
    requireNamespace("ggplot2")
    
    if (is.null(obs)) {
        obs <- scKirby::get_obs(obj = obj,
                                verbose = verbose)
    }
    if (is.null(rownames(obs))) {
        stopper("obs must have rownames to merge with obj.")
    } 
    #### Prepare data #### 
    obsm <- scKirby::get_obsm(obj = obj,
                              keys = keys) 
    keys <- names(obsm) 
    obsm <- scKirby::get_n_elements(l = obsm,
                                    n = 1,
                                    verbose = verbose) 
    if(x_dim>ncol(obsm) || x_dim<1) {
        stopper(
            "x_dim must be an integer >= the number of obsm dimensions."
            ) 
    }
    if(y_dim>ncol(obsm) || y_dim<1) {
        stopper(
            "y_dim must be an integer >= the number of obsm dimensions."
            ) 
    }
    #### Fix rownames to match obs rownames ####
    if (isTRUE(fix_rownames)) {
        obs <- obs |>
            `rownames<-`(gsub("-|[:]|[.]", "_", rownames(obs)))
        obsm <- obsm |>
            `rownames<-`(gsub("-|[:]|[.]", "_", rownames(obsm)))
    }
    keep_cols <- grep(paste(paste0("^",c("umap",keys)),collapse = "|"),
                      names(obs),ignore.case = TRUE, invert = TRUE)
    plot_df <- merge(
        x = obsm,
        y = obs[,keep_cols],
        by = 0
    )
    if(nrow(plot_df)==0){
        stopper("No rows in merged obs + obsm data.")
    } else {
        messager("Constructed plot_df:",
                 formatC(nrow(plot_df),big.mark = ","),"rows x",
                 formatC(ncol(plot_df),big.mark = ","),"columns"
                 )
    }
    #### Plot ####
    gp <- ggplot2::ggplot(plot_df, ggplot2::aes_string(
        x = colnames(obsm)[x_dim],
        y = colnames(obsm)[y_dim],
        color = color_var, 
        label = label_var
    )) +
        ggplot2::geom_point(ggplot2::aes_string(size = size_var),
                            alpha = point_alpha) +
        ggplot2::theme_bw()
    if (isTRUE(labels)) {
        gp <- gp + ggrepel::geom_text_repel(max.overlaps = max.overlaps)
    }
    if (isTRUE(show_plot)) methods::show(gp)
    return(gp)
}
