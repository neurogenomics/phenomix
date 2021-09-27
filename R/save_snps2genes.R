save_snps2genes <- function(gene_hits,
                            save_dir,
                            sumstats_file,
                            method = NULL,
                            nCores = 1,
                            verbose = TRUE) {
    if (!is.null(save_dir)) {
        #### Construct path ####
        if (methods::is(sumstats_file, "character")) {
            id <- gsub(
                paste(MungeSumstats:::supported_suffixes(), collapse = "|"),
                "", basename(sumstats_file)
            )
            out_path <- file.path(save_dir, paste0(id, "_", method, ".tsv.gz"))
        } else {
            out_path <- paste0(tempfile(), method, ".tsv.gz")
        }
        #### Save ####
        messager("Saving gene_hits ==>", out_path, v = verbose)
        data.table::fwrite(gene_hits, out_path, nThread = nCores)
        return(out_path)
    } else {
        return(NULL)
    }
}
