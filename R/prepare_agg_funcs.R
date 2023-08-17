#' prepare_agg_funcs
#' 
#' @keywords internal
#' @importFrom stats sd
#' @importFrom data.table uniqueN
prepare_agg_funcs <- function(merged_hits,
                              extra_cols = c("ABC.Score", "CellType"),
                              verbose = TRUE) {
    #### Define functions ####
    mean_ <- function(x) {
        mean(as.numeric(x), na.rm = TRUE)
    }
    sd_ <- function(x) {
        stats::sd(as.numeric(x), na.rm = TRUE)
    }
    uniqueN_ <- function(x) {
        data.table::uniqueN(x, na.rm = TRUE)
    }
    unique_ <- function(x) {
        unique(x)[1]
    }
    N_ <- function(x) {
        length(x)
    }
    collapse_ <- function(x) {
        paste(unique(x), collapse = ";")
    }
    #### Define columns ####
    agg_functions <- list(
        CHR = unique_, BP_min = min,BP_max = max,
        BETA_mean = mean_, P_mean = mean_, SE_mean = mean_,
        BETA_sd = sd_, P_sd = sd_, SE_sd = sd_,
        ROWS = N_, NSNPS = uniqueN_
    )
    col_selection <- c(
        "CHR","BP","BP",
        "BETA", "P", "SE",
        "BETA", "P", "SE",
        "SNP", "SNP"
    )
    i <- which(col_selection %in% names(merged_hits))
    agg_functions <- agg_functions[i]
    col_selection <- col_selection[i]
    #### Add extra columns ####
    extra_cols <- extra_cols[extra_cols %in% names(merged_hits)]
    if (length(extra_cols) > 0) {
        messager("Extra cols used:", paste(extra_cols, collapse = ", "),
            v = verbose
        ) 
        if ("ABC.Score" %in% extra_cols) {
            agg_functions[[paste0("ABC.Score", "_mean")]] <- mean_
            col_selection <- c(col_selection, "ABC.Score")
            agg_functions[[paste0("ABC.Score", "_sd")]] <- sd_
            col_selection <- c(col_selection, "ABC.Score")
        }
        if ("CellType" %in% extra_cols) {
            agg_functions[["CellTypes"]] <- collapse_
            col_selection <- c(col_selection, "CellType")
        }
    }
    return(list(
        agg_functions = agg_functions,
        col_selection = col_selection
    ))
}
