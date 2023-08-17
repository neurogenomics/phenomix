get_gene_length <- function(gene_hits,
                            ref_genome = "GRCh37",
                            gene_col = "SYMBOL",
                            use_symbols = TRUE,
                            drop_na = TRUE,
                            verbose = TRUE) {
    # devoptera::args2vars(get_gene_length) 
    
    messager("Computing gene lengths.", v = verbose)
    if(!methods::is(gene_hits,"data.table")){
        gene_hits <- data.table::data.table(gene_hits)
    } 
    txdb <- select_txdb_build(
        ref_genome = ref_genome,
        verbose = verbose
    )
    #### Get transcript ranges ####
    length_key <- GenomicFeatures::genes(txdb,
        columns = c(
            "GENEID", "TXID",
            "TXSTART", "TXEND"
        )
    ) 
    length_key$GENEID <- as.character(length_key$GENEID)
    length_key <- data.table::data.table(data.frame(length_key),
                                         key = "GENEID")
    #### Translate gene IDs back to gene symbols ####
    if (gene_col!="GENEID" &&
        isTRUE(use_symbols)) {
        length_key <- translate_geneids_txdb(
            gene_hits = length_key,
            gene_col = "GENEID",
            verbose = verbose
        )
        gene_col <- "SYMBOL"
        data.table::setkeyv(length_key, "SYMBOL")
    } else {
        data.table::setkeyv(length_key, "GENEID")
    }
    #### Query length_key for GENELEN ####
    gene_hits <- merge(
        gene_hits,
        length_key[,c(gene_col,"start","end","width"), with=FALSE] |>
            data.table::setnames(c("start","end","width"),
                                 c("GENE_start","GENE_end","GENELEN"),
                                 skip_absent = TRUE),
        by = gene_col, 
        all.x = TRUE, 
        sort = FALSE
    ) 
    #### Drop NAs ####
    na_count <- sum(is.na(gene_hits$GENELEN))
    if (isTRUE(drop_na) &&
        na_count > 0) {
        messager("Dropping", formatC(na_count, big.mark = ","),
            "genes without GENELEN.",
            v = verbose
        )
        gene_hits <- stats::na.omit(gene_hits, cols = "GENELEN")
    }
    return(gene_hits)
}
