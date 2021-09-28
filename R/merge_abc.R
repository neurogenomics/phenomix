merge_abc <- function(sumstats,
                      dataset = "Nasser2020",
                      abc = NULL,
                      nCores = 1,
                      verbose = TRUE) {
    CHR <- NULL;
    #### Import ABC model predictions ####
    abc <- import_abc(
        dataset = dataset,
        abc = abc,
        nCores = nCores,
        verbose = verbose
    )
    #### data.table method ####
    if (!chr_prefix(sumstats)) {
        sumstats[, CHR := paste0("chr", CHR)]
    }
    data.table::setkeyv(sumstats, c("CHR", "BP"))
    ##### Merge ####
    messager("Merging sumstats with ABC model predictions.", v = verbose)
    sumstats_abc <- merge(x = sumstats, y = abc, all.x = FALSE, all.y = FALSE)
    #### Sort ####
    messager("Sorting genomic coordinates.", v = verbose)
    sumstats_abc[, CHR := as.integer(gsub("chr", "", CHR))]
    data.table::setkeyv(sumstats_abc, c("CHR", "BP"))
    #### Report ####
    messager(paste0(
        "Returning merged sumstats-ABC data.table with:\n",
        paste0(" - ", formatC(data.table::uniqueN(sumstats_abc$CellType), big.mark = ","), " chromosomes."), "\n",
        paste0(" - ", formatC(nrow(sumstats_abc), big.mark = ","), " rows."), "\n",
        paste0(" - ", formatC(data.table::uniqueN(sumstats_abc$SNP), big.mark = ","), " SNPs."), "\n",
        paste0(" - ", formatC(data.table::uniqueN(sumstats_abc$CellType), big.mark = ","), " CellTypes."), "\n",
        paste0(" - ", formatC(data.table::uniqueN(sumstats_abc$TargetGene), big.mark = ","), " TargetGenes."), "\n",
        paste0(" - ", formatC(data.table::uniqueN(sumstats_abc$class), big.mark = ","), " genomic region classes.")
    ), v = verbose)
    #### Return ####
    return(sumstats_abc)
}
