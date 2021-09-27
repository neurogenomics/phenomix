#' Translate GENEId to SYMBOL
#'
#' Translate \code{GENEID} from
#' \link[TxDb.Hsapiens.UCSC.hg19.knownGene]{TxDb.Hsapiens.UCSC.hg19.knownGene}
#' or
#' \link[TxDb.Hsapiens.UCSC.hg38.knownGene]{TxDb.Hsapiens.UCSC.hg38.knownGene}
#' to gene \code{SYMBOL}.
#'
#' @param gene_hits Output from \code{aggregate_sumstats}.
#' @param gene_var Column name with gene IDs.
#' @param drop_na Drop genes without a corresponding \code{SYMBOL}
#' @param verbose Print messages.
#'
#' @keywords internal
#' @importFrom AnnotationDbi select
#' @importFrom data.table data.table setkey
#' @importFrom dplyr %>%
#' @importFrom stats na.omit
translate_geneids_txdb <- function(gene_hits,
                                   gene_var = "GENEID",
                                   drop_na = TRUE,
                                   verbose = TRUE) {
    messager("Translating", gene_var, "to SYMBOL.", v = verbose)
    symbol_key <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
        keys = gene_hits[[gene_var]],
        columns = c("SYMBOL", "ENTREZID"),
        keytype = "ENTREZID"
    ) %>%
        data.table::data.table(key = "ENTREZID") %>%
        na.omit() %>%
        unique()
    data.table::setkeyv(gene_hits, gene_var)
    gene_hits[, SYMBOL := symbol_key[get(gene_var), "SYMBOL"]]
    #### Drop NAs ####
    na_count <- sum(is.na(gene_hits$SYMBOL))
    if (drop_na & na_count > 0) {
        messager("Dropping", na_count, "genes without SYMBOL.",
            v = verbose
        )
        gene_hits <- stats::na.omit(gene_hits, "SYMBOL")
    }
    return(gene_hits)
}
