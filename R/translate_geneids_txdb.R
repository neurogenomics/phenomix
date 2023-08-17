#' Translate GENEId to SYMBOL
#'
#' Translate \code{GENEID} from
#' \link[TxDb.Hsapiens.UCSC.hg19.knownGene]{TxDb.Hsapiens.UCSC.hg19.knownGene}
#' or
#' \link[TxDb.Hsapiens.UCSC.hg38.knownGene]{TxDb.Hsapiens.UCSC.hg38.knownGene}
#' to gene \code{SYMBOL}.
#'
#' @param gene_hits Output from \code{aggregate_sumstats}.
#' @param gene_col Column name with gene IDs.
#' @param drop_na Drop genes without a corresponding \code{SYMBOL}
#' @param verbose Print messages.
#'
#' @keywords internal
#' @importFrom AnnotationDbi select
#' @importFrom data.table data.table setkey 
#' @importFrom stats na.omit
translate_geneids_txdb <- function(gene_hits,
                                   gene_col = "GENEID",
                                   drop_na = TRUE,
                                   verbose = TRUE) {
    SYMBOL <- NULL;
    messager("Translating", gene_col, "to SYMBOL.", v = verbose)
    symbol_key <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
        keys = as.character(gene_hits[[gene_col]]),
        columns = c("SYMBOL", "ENTREZID"),
        keytype = "ENTREZID"
    ) |>
        data.table::data.table(key = "ENTREZID") |>
        na.omit() |>
        unique()
    data.table::setkeyv(gene_hits, gene_col)
    #### Make sure input/output cols are unique ####
    if(gene_col=="SYMBOL"){
        messager("Renaming gene_col as SYMBOL.1",v=verbose)
        data.table::setnames(gene_hits,gene_col,"SYMBOL.1")
    }
    gene_hits[, SYMBOL := symbol_key[get(gene_col), "SYMBOL"]]
    #### Drop NAs ####
    na_count <- sum(is.na(gene_hits$SYMBOL))
    if (isTRUE(drop_na) && 
        na_count > 0) {
        messager("Dropping",
                 formatC(na_count,big.mark = ","),
                 "rows without SYMBOL.",
            v = verbose
        )
        gene_hits <- stats::na.omit(gene_hits, "SYMBOL")
    } 
    return(gene_hits)
}
