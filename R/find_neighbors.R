#' Search nearest neighbors
#'
#' Search [shared] K-nearest neighbor graph to find the
#'  samples that are most similar to those matching a substring search.
#'
#' @param obj \pkg{Seurat} object.
#' @param assay Assay to use.
#' @param slot Slot to use. 
#' @param graph_key Name of the graph to use.
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
#' @inheritParams scKirby::get_obsm
#'
#' @export 
#' @importFrom stats setNames  
#' @import data.table 
#' @examples 
#' obj <- get_HPO()[seq(100),]
#' top_neighbors <- find_neighbors(obj = obj,
#'                                 var1_search = "parkinson",
#'                                 label_col = "HPO_label")
find_neighbors <- function(obj,
                           graph_key = NULL,
                           assay = NULL,
                           slot = NULL,
                           keys = NULL,
                           var1_search = NULL,
                           label_col = NULL,
                           var2_group = NULL,
                           group_col = NULL,
                           max_neighbors = Inf,
                           add_original_names = TRUE,
                           verbose = TRUE) {
    
    # devoptera::args2vars(find_neighbors); slot=NULL;
    
    trait1 <- trait2 <- similarity <- NULL;
    #### Get obs ####
    obs <- scKirby::get_obs(obj = obj,
                            verbose = verbose)
    if(!is.null(label_col)){
        if(!label_col %in% names(obs)){
            stopper(paste0("label_col=",shQuote(label_col)),
                    "not found in obs metadata.")
        }
    }
    #### Extract/compute correlation matrix ####
    cmat <- get_cor(
        obj = obj,
        graph_key = graph_key, 
        keys = keys,
        assay = assay,
        slot = slot,
        return_obj = FALSE,
        verbose = verbose
    ) 
    #### Align naming of obs metadata and matrix ####
    if (is.null(label_col)) {
        sample_names <- rownames(cmat)
    } else {
        sample_names <- make.unique(obs[[label_col]])
    }
    obs$sample_names <- sample_names
    names_dict <- stats::setNames(rownames(cmat), sample_names) 
    fgraph <- cmat |>
        `row.names<-`(sample_names) |>
        `colnames<-`(sample_names)
    if(is.null(var1_search)){
        stopper("Must provide terms to var1_search.")
    } else {
        messager("Filtering results by var1_search:",
            paste(var1_search,collapse = " | "),
            v = verbose
        )
        targets1 <- grep(
            pattern = paste(var1_search,collapse="|"),
            x = sample_names,
            ignore.case = TRUE, value = TRUE
        )
        if (length(targets1) > 0) {
            messager("+", formatC(length(targets1), big.mark = ","),
                "entries matching var1_search identified.",
                v = verbose
            ) 
                #### Ensure the rest doesn't break when only 1 var1 hit ####
                # Matrices get converted to vectors if there's only 1 col  
            if(length(targets1)!=1){
                fgraph <- fgraph[targets1,,drop=FALSE]
            } 
        } else {
            stopper(
                "0 entries in label_col match the ",
                "substring search for `var1_search`"
            )
        }
    }

    if (!is.null(group_col) &&
        (group_col %in% colnames(obs)) &&
        (!is.null(var2_group))) {
        messager("+ Filtering results by var2_group:", var2_group,
                 v = verbose)
        targets2 <- obs[
            grepl(pattern = paste(var2_group, collapse = "|"), 
                  x = obs[[group_col]],
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
            stopper(
                "0 entries in group_col match the ",
                "substring search for var2_group."
            )
        }
    } 

    
    messager("Converting graph to data.table",v=verbose) 
    top_candidates <- fgraph[targets1,,drop=FALSE] |>
        data.table::as.data.table(keep.rownames="trait1") |>
        data.table::melt.data.table(id.vars="trait1",
                                    variable.name="trait2",
                                    value.name = "similarity", 
                                    drop=FALSE)
    top_candidates[,trait1:=as.character(trait1)][,trait2:=as.character(trait2)]
    top_candidates <- top_candidates[trait1 != trait2][similarity>0] 
    top_candidates <- 
        top_candidates[,.SD[similarity %in% head(sort(similarity),
                                                 max_neighbors)],
                       by="trait1"] |>
        data.table::setorderv("similarity",-1)
         
    
    if (!is.null(group_col)) {
        messager("+ Adding group_col to results.", v = verbose)
        keys <- stats::setNames(obs[[group_col]], obs$sample_names)
        top_candidates$trait2_group <- keys[top_candidates$trait2]
    }
    if (isTRUE(add_original_names)) {
        messager("+ Adding original names to results.", v = verbose)
        top_candidates$trait1_id <- names_dict[top_candidates$trait1]
        top_candidates$trait2_id <- names_dict[top_candidates$trait2]
    }
    messager("+ Returning", formatC(nrow(top_candidates), big.mark = ","),
             "pair-wise similarities.",v = verbose)
    return(top_candidates)
}
