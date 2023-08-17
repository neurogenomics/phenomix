remove_mhc_genes <- function(dat,
                             gene_col = "SYMBOL",
                             verbose = TRUE,
                             ...) { 
    MHC_genes <- get_mhc_genes(
        verbose = verbose,
        ...
    ) 
    #### Check which MHC data column to filter by  #### 
    ## Entrez gene ID or HGNC gene symbol
    MHC_col <-if(methods::is(dat[[gene_col]],"numeric")){
        "GENEID"
    } else {
        "SYMBOL"
    }  
    #### Filter
    MHC_hits <- unique(
        dat[get(gene_col) %in% unique(MHC_genes[[MHC_col]]),][[gene_col]]
    )  
    messager("Removing",formatC(length(MHC_hits), big.mark = ","), 
             "MHC gene(s) from data.",v = verbose)
    dat <- dat[!get(gene_col) %in% MHC_hits, ]
    return(dat)
}
