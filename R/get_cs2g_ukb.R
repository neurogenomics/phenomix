#' Import precomputed cS2G results
#'
#' Import 47 GWAS (from UK Biobank) mapped to genes using the
#' cS2G strategy.
#'
#' @section References:
#'
#' \href{https://doi.org/10.1101/2021.08.02.21261488}{
#' Gazal et al. 2021 preprint}
#'
#' \href{https://alkesgroup.broadinstitute.org/cS2G}{
#' Broad Institute server}
#'
#' @param URL URL to cS2G dataset.
#' @param as_matrix Convert the data to a sparse matrix.
#' @param nThread Number of threads to use when reading in data.
#' @param verbose Print messages.
#'
#' @return Named with two items: "data" and "metadata"
#'
#' @export
#' @importFrom data.table fread dcast
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr %>%
#' @import Matrix
get_cs2g_ukb <- function(URL = file.path(
                             "https://storage.googleapis.com",
                             "broad-alkesgroup-public/cS2G",
                             "finemapping_cS2G_UKBB",
                             "finemapping_cS2G_UKBB_PIP05.annot"
                         ),
                         as_matrix = TRUE,
                         value.var = c("cS2G", "PIP"),
                         nThread = 1,
                         verbose = TRUE) {
    value.var <- value.var[1]
    dat <- data.table::fread(URL, nThread = nThread)
    #### Prepare metadata ####
    dat[,ID:=gsub(" ","-",DISEASE.TRAIT)]
    #### Remove non-ASCI characters ####
    dat[, ID:=make.unique(iconv(ID, "latin1", "ASCII", sub=""),sep = "_")]
   
    metadata <- dat[, .(
        DISEASE.TRAIT = unique(DISEASE.TRAIT),
        N_SNP = data.table::uniqueN(SNP),
        N_GENE = data.table::uniqueN(gene),
        cS2G_mean = mean(cS2G = TRUE),
        PIP_mean = mean(PIP, na.rm = TRUE)
    ),
    by = c("ID")
    ]
    metadata[, GROUP := stringr::str_split(DISEASE.TRAIT, "_", 
                                           n = 2, simplify = TRUE)[, 1]]
    metadata[, TRAIT := stringr::str_split(DISEASE.TRAIT, "_", 
                                           n = 2, simplify = TRUE)[, 2]]
    metadata <- data.frame(metadata, row.names = metadata$ID)


    if (as_matrix) {
        messager("Converting to matrix format.", v = verbose)
        mat <- data.table::dcast(
            data = dat,
            formula = gene ~ ID,
            value.var = value.var,
            fill = 0,
            fun.aggregate = "mean"
        ) %>%
            tibble::column_to_rownames("gene") %>%
            as.matrix() %>%
            methods::as("sparseMatrix")
        #### Ensure rownames match ####
        metadata <- metadata[colnames(mat), ]
        # sum(rownames(metadata)!=colnames(mat))
        return(list(
            data = mat,
            metadata = metadata
        ))
    } else {
        return(list(
            data = dat,
            metadata = metadata
        ))
    }
}
