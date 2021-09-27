#' Plot \pkg{variancePartition}
#'
#' Plot \link[variancePartition]{fitExtractVarPartModel}
#' output.
#'
#' @param varPart \link[variancePartition]{fitExtractVarPartModel}
#' output.
#' @param plot_PercentBars Add  \link[variancePartition]{plotPercentBars}.
#' @param plot_VarPart Add  \link[variancePartition]{plotPercentBars}.
#' @param genes Genes to plot in \link[variancePartition]{plotVarPart}.
#' @param n_genes Number of genes to select when \code{genes=NULL}.
#' @param ... Additional arguments passed to
#'  \link[patchwork]{wrap_plots}.
#'
#' @importFrom variancePartition sortCols plotPercentBars  plotVarPart
#'
#' @export
plot_variancePartition <- function(varPart,
                                   plot_PercentBars = TRUE,
                                   plot_VarPart = TRUE,
                                   genes = NULL,
                                   n_genes = 10,
                                   show_plot = TRUE,
                                   ...) {
    out <- list()
    vp <- variancePartition::sortCols(varPart)
    #### Select genes ####
    if (is.null(genes)) {
        message(
            "Selecting ", n_genes,
            " with the greatest total variance explained by metadata variables."
        )
        genes <- names(tail(sort(rowSums(vp[, -ncol(vp)])), n_genes))
    } else {
        genes <- genes[genes %in% rownames(vp)]
    }
    #### plotPercentBars ####
    if (plot_PercentBars) {
        ppb <- variancePartition::plotPercentBars(vp[genes, ])
        out[["plotPercentBars"]] <- ppb
    }
    #### plotVarPart ####
    if (plot_VarPart) {
        vp2 <- variancePartition::sortCols(varPart, decreasing = FALSE)
        pvp <- variancePartition::plotVarPart(vp2,
            label.angle = 45
        ) +
            stat_summary(aes(label = ..y..),
                fun = function(x) {
                    round(median(x), 2)
                },
                geom = "label", size = 3, alpha = .75
            )
        out[["plotVarPart"]] <- pvp
    }

    out_merged <- patchwork::wrap_plots(out, ...)
    if (show_plot) print(out_merged)
    return(out_merged)
}
