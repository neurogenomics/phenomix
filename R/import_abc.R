#' Import ABC model predictions
#'
#' Import Activity-By-Contact (ABC) model predictions
#' from multiple publications.
#'
#' @section References:
#'
#' \bold{Nasser2020}\cr
#' \href{https://www.nature.com/articles/s41586-021-03446-x}{Nasser et al. 2020, Nature}
#'
#' \bold{Fulco2019}\cr
#' \href{https://www.nature.com/articles/s41588-019-0538-0?proof=t}{Fulco et al. 2019, Nature}
#' \href{https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction}{Fulco et al. 2019, GitHub}
#'
#' \href{https://www.engreitzlab.org/resources/}{ABC models from multiple Engreitz Lab publications}
#'
#' @param dataset Which ABC model to import.
#' @param nCores Number of cores to parallelise across. 
#' @param abc Use a previously downloaded ABC \link[data.table]{data.table}.
#' @param save_dir Directory to cache files to.
#' @param verbose Print messages.
#'
#' @return \code{gene_hits} \link[data.table]{data.table}
#'
#' @export
#' @importFrom data.table setkey fread
import_abc <- function(dataset = "Nasser2020",
                       abc = NULL,
                       nCores = 1,
                       save_dir = tools::R_user_dir(
                           package = "phenomix",
                           which = "cache"),
                       verbose = TRUE) {
    dataset <- tolower(dataset[1])
    
    #### Return provided data ####
    if (!is.null(abc)) {
        messager("Using supplied ABC data.table", v = verbose)
        return(abc)
    }
    #### NOTE: aligned to HG19 ####
    if (dataset == "nasser2020") {
        abc_url <- paste0(
            "ftp://ftp.broadinstitute.org/outgoing/",
            "lincRNA/ABC/",
            "AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz"
        )
    }
    #### Read ABC ####
    path <- file.path(save_dir,basename(abc_url))
    if(!file.exists(path)){
        dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
        utils::download.file(url = abc_url,
                             destfile = path)
    }
    messager("Importing",dataset,"ABC model predictions.",
        v = verbose
    )
    abc <- data.table::fread(path, nThread = nCores)
    data.table::setnames(abc, c("chr", "start"), c("CHR", "BP"))
    data.table::setkeyv(abc, c("CHR", "BP"))
    return(abc)
}
