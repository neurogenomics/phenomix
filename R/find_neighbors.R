#' Search nearest neighbors
#'
#' Search [shared] K-nearest neighbor graph to find the
#'  samples that are most similar to those matching a substring search.
#'
#' @param obj \pkg{Seurat} object.
#' @param assay Assay to use.
#' @param slot Slot to use.
#' @param reduction Reduction to use.
#' @param graph_name Name of the graph to use.
#' If none provided, will use the last graph available.
#' If no graphs are available, new ones will be computed using
#'  Pearson's correlation on the raw counts matrix.
#' @param var1_search Substring search term to filter var1 by.
#' If a vector is supplied instead, this will be interpreted as an "or" query.
#' @param label_col \code{meta.data} column used to name the rows/columns of the graph.
#' \code{label_col} will also be used in the search for \code{var1_search} substring.
#' @param var2_group Substring search term to filter var2 by,
#' according to  \code{group_col}.
#' @param group_col \code{meta.data} column used to filter var2
#' when \code{var2_group} is used.
#' @param max_neighbors The max number of neighbors (var2) per term (var1).
#' @param add_original_names Add original names into the results.
#' This can be useful when var1 names are forced to be unique internally.
#' @param verbose Whether to print messages.
#'
#' @examples 
#' degas <- get_DEGAS()
#'
#' ### No group filter
#' top_neighbors <- find_neighbors(
#'     obj = degas,
#'     var1_search = "parkinson", 
#'     label_col = "label_phe"
#' )
#' @export
#' @importFrom dplyr %>% mutate_at arrange group_by slice_max desc rename
#' @importFrom stats setNames
#' @importFrom Matrix as.matrix
#' @importFrom reshape2 melt
#' @importFrom data.table data.table
find_neighbors <- function(obj,
                           graph_name = NULL,
                           assay = NULL,
                           slot = NULL,
                           reduction=NULL,
                           var1_search = NULL,
                           label_col = NULL,
                           var2_group = NULL,
                           group_col = NULL,
                           max_neighbors = Inf,
                           add_original_names = TRUE,
                           verbose = TRUE) {
    Var1 <- Var2 <- trait1 <- trait2 <- similarity <- NULL;
    #### Extract/compute correlation matrix ####
    cmat <- extract_cor(
        obj = obj,
        graph_name = graph_name,
        assay = assay,
        slot = slot,
        reduction = reduction,
        verbose = verbose
    )
    metadata <- extract_metadata(obj = obj)
    if (is.null(label_col)) {
        sample_names <- rownames(cmat)
    } else {
        sample_names <- make.unique(metadata[[label_col]])
    }
    metadata$sample_names <- sample_names
    names_dict <- stats::setNames(rownames(cmat), sample_names)

    fgraph <- cmat %>%
        `row.names<-`(sample_names) %>%
        `colnames<-`(sample_names)

    if (!is.null(var1_search)) {
        messager("+ Filtering results by var1_search:",
            var1_search,
            v = verbose
        )
        targets1 <- grep(
            pattern = var1_search, x = sample_names,
            ignore.case = TRUE, value = TRUE
        )
        if (length(targets1) > 0) {
            messager("+", formatC(length(targets1), big.mark = ","),
                "entries matching var1_search identified.",
                v = verbose
            )
            fgraph <- fgraph[targets1, ]
        } else {
            stop(
                "0 entries in label_col match the ",
                "substring search for `var1_search`"
            )
        }
    }

    if (!is.null(group_col) && (group_col %in% colnames(metadata))) {
        messager("+ Filtering results by var2_group:", var2_group, v = verbose)
        targets2 <- metadata[
            grepl(var2_group, metadata[[group_col]],
                ignore.case = TRUE
            ),
            "sample_names"
        ]
        if (length(targets2) > 0) {
            messager("+", formatC(length(targets2), big.mark = ","),
                "entries matching var2_group identified.",
                v = verbose
            )
            fgraph <- fgraph[, targets2]
        } else {
            stop(
                "0 entries in group_col match the ",
                "substring search for var2_group."
            )
        }
    }

    top_candidates <- fgraph %>%
        Matrix::as.matrix() %>%
        reshape2::melt(value.name = "similarity") %>%
        dplyr::rename(trait1 = Var1, trait2 = Var2) %>%
        data.table::data.table() %>%
        dplyr::mutate_at(c("trait1", "trait2"), as.character) %>%
        subset(trait1 != trait2) %>%
        subset(similarity > 0) %>%
        dplyr::group_by(trait1) %>%
        dplyr::slice_max(
            order_by = similarity,
            n = max_neighbors
        ) %>%
        dplyr::arrange(dplyr::desc(similarity)) %>%
        data.table::data.table()

    if (!is.null(group_col)) {
        messager("+ Adding group_col to results.", v = verbose)
        keys <- stats::setNames(metadata[[group_col]], metadata$sample_names)
        top_candidates$trait2_group <- keys[top_candidates$trait2]
    }
    if (add_original_names) {
        messager("+ Adding original names to results.", v = verbose)
        top_candidates$trait1_id <- names_dict[top_candidates$trait1]
        top_candidates$trait2_id <- names_dict[top_candidates$trait2]
    }
    messager("+ Returning", formatC(nrow(top_candidates), big.mark = ","),
        "pair-wise similarities.",
        v = verbose
    )
    return(top_candidates)
}
