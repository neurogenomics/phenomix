save_snps2genes <- function(gene_hits,
                            save_path,
                            nCores = 1,
                            verbose = TRUE) {
    if (!is.null(save_path)) { 
        #### Save ####
        messager("Saving gene_hits ==>", save_path, v = verbose)
        dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
        data.table::fwrite(gene_hits, save_path, nThread = nCores)
        return(save_path)
    } else {
        return(NULL)
    }
}
