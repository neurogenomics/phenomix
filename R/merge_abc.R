merge_abc <- function(dat,
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
    if (!chr_prefix(dat)) {
        dat[, CHR:=paste0("chr",CHR)]
    }
    data.table::setkeyv(dat, c("CHR", "BP"))
    ##### Merge ####
    messager("Merging dat with ABC model predictions.", v = verbose)
    dat_abc <- merge(x = dat, 
                     y = abc,
                     all.x = FALSE,
                     all.y = FALSE, 
                     by = c("CHR","BP"))
    if(nrow(dat_abc)==0){
        stopper("No rows in dat could be annotated with ABC dataset.")
    }
    #### Sort ####
    messager("Sorting genomic coordinates.", v = verbose)
    dat_abc[, CHR:=as.integer(gsub("chr","",CHR))]
    data.table::setkeyv(dat_abc, c("CHR", "BP"))
    #### Report ####
    messager(paste0(
        "Returning merged dat-ABC data.table with:\n",
        paste0(" - ", formatC(data.table::uniqueN(dat_abc$CellType), big.mark = ","), " chromosomes."), "\n",
        paste0(" - ", formatC(nrow(dat_abc), big.mark = ","), " rows."), "\n",
        paste0(" - ", formatC(data.table::uniqueN(dat_abc$SNP), big.mark = ","), " SNPs."), "\n",
        paste0(" - ", formatC(data.table::uniqueN(dat_abc$CellType), big.mark = ","), " CellTypes."), "\n",
        paste0(" - ", formatC(data.table::uniqueN(dat_abc$TargetGene), big.mark = ","), " TargetGenes."), "\n",
        paste0(" - ", formatC(data.table::uniqueN(dat_abc$class), big.mark = ","), " genomic region classes.")
    ), v = verbose)
    #### Return ####
    return(dat_abc)
}
