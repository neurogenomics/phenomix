#' Import precomputed cS2G results
#'
#' Import 5,981 GWAS (from GWAS Catalog) mapped to genes using the
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
#' @importFrom data.table fread dcast .GRP
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr %>% 
get_cs2g_gwascatalog <- function(URL = file.path(
                                     "https://storage.googleapis.com",
                                     "broad-alkesgroup-public/cS2G",
                                     "gwas_catalog_cS2G",
                                     "gwas_catalog_v1.0-associations_e100_r2021-02-25.annot"
                                 ),
                                 as_matrix = TRUE,
                                 nThread = 1,
                                 verbose = TRUE) {
    DISEASE.TRAIT <- . <- PUBMEDID <- ID <- SNP <- gene <- 
        cS2G <- trait_id <- NULL;
    dat <- data.table::fread(URL, nThread = nThread)
    #### Assign unique trait IDs to differentiate them after truncation ####
    dat[, trait_id := .GRP, by = .(DISEASE.TRAIT, PUBMEDID)]
    #### Remove non-ASCI characters ####
    dat[, ID :=
        paste(
            stringr::str_trunc(
                iconv(
                    paste(
                        gsub(" ", "-", DISEASE.TRAIT),
                        PUBMEDID,
                        sep = "_"
                    ),
                    "latin1", "ASCII",
                    sub = ""
                ),
                width = 40
            ), # MOFA2 doesn't let you have IDs > 50 characters
            trait_id,
            sep = "_"
        )]


    metadata <- dat[, .(
        N_SNP = data.table::uniqueN(SNP),
        N_GENE = data.table::uniqueN(gene),
        cS2G_mean = mean(cS2G)
    ),
    by = c("ID", "DISEASE.TRAIT", "PUBMEDID")
    ]
    metadata <- data.frame(metadata, row.names = metadata$ID)

    if (as_matrix) {
        messager("Converting to matrix format.", v = verbose)
        mat <- data.table::dcast(
            data = dat,
            formula = gene ~ ID,
            value.var = "cS2G",
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
