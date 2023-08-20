#' Map SNPs to genes
#'
#' Map SNPs (both coding and non-coding) onto genes
#' using different methods.
#'
#' @section Other packages worth exploring:
#' \href{https://github.com/thewonlab/H-MAGMA}{H-MAGMA}
#' Really hard to use with no docs on how to create new annot files.
#' It's literally just MAGMA (when you provide some extra args).
#'
#' \href{https://www.rdocumentation.org/packages/PAGWAS/versions/2.0/topics/snps.to.genes}{PAGWAS}
#'
#' \href{https://cran.r-project.org/web/packages/rsnps/rsnps.pdf}{rsnp}
#'
#' \href{orthogene}{https://github.com/neurogenomics/orthogene}
#' Protein-coding variants only?
#'
#' \href{OpenTargets}{https://genetics.opentargets.org/}
#' Could use colocalization results, but these are both GWAS- and QTL-specific.
#' Could compute mean coloc score per gene for each SNP and use that.
#' Other position-based functional data might still be useful.
#'
#' @param dat GWAS summary statistics munged by
#' \link[MungeSumstats]{format_sumstats}.
#' Can be a path to the saved file or \link[data.table]{data.table}.
#' @param adjust_z Whether to adjust Z-statistic using
#'  the \code{model} and \code{formula}.
#' @param nCores Number of cores to parallelise across.
#' @param method Method to use for mapping SNPs to genes.
#' @param save_path Path to save results to.
#' @param verbose Print messages.
#' @param return_path Whether to return the path to the saved results.
#' @inheritParams import_abc
#' @inheritParams adjust_zstat
#' @inheritParams aggregate_sumstats
#' @inheritParams translate_geneids_txdb
#'
#' @return gene_hits \link[data.table]{data.table}.
#'
#' @export
#' @import data.table
#' @examples
#' dat <- MungeSumstats::formatted_example()
#' dat2 <- map_snps2genes(dat)
map_snps2genes <- function(dat,
                           #### Methods args ####
                           method = c("txdb","abc"),
                           dataset = "Nasser2020",
                           abc = NULL,
                           agg_var = "SYMBOL",
                           drop_na = TRUE,
                           #### Adjust Z-stat args ####
                           adjust_z = TRUE,
                           drop_mhc = TRUE,
                           model = NULL,
                           log_vars = c("NSNPS", "NPARAM", "GENELEN"),
                           formula = ZSTAT ~ NSNPS + logNSNPS + NPARAM +
                               logNPARAM + GENELEN + logGENELEN,
                           #### Util args ####
                           nCores = 1,
                           save_path = tempfile(),
                           return_path = FALSE,
                           verbose = TRUE) {
    # devoptera::args2vars(map_snps2genes)
    
    check_required_cols(dat = dat)
    start <- Sys.time()
    #### Check method ####
    method <- tolower(method[1])
    #### Prepare dat ####
    # Replace with MungeSumstats:::format_sumstats() once it's fixed.
    if (methods::is(dat, "character")) {
        messager("Importing dat.", v = verbose)
        dat <- data.table::fread(dat,
                                 nThread = nCores)
    } else {
        dat <- data.table::as.data.table(dat)
    }
    #### Select method and map SNPs ####
    #### txdb ####
    if (method == "txdb") {
        #### Map #### 
        gene_hits <- map_snps2genes_txdb(
            dat = dat,
            agg_var = agg_var,
            drop_na = drop_na,
            nCores = nCores,
            verbose = verbose
        ) 
    }
    #### ABC ####
    if (method == "abc") {
        #### Map #### 
        gene_hits <- map_snps2genes_abc(
            dat = dat,
            abc = abc,
            agg_var = agg_var,
            dataset = dataset,
            nCores = nCores,
            verbose = verbose
        ) 
    }
    #### Adjust ZSTAT ####
    if (isTRUE(adjust_z)) {
        gene_hits <- adjust_zstat(
            dat = gene_hits,
            model = model,
            log_vars = log_vars,
            drop_mhc = drop_mhc,
            formula = formula,
            verbose = verbose
        ) 
    }
    #### OpenTargets ####
    # if(method=="opentargets"){
    #     # Query OT for coloc probabilities for each locus
    #     # and then convert to genome-wide gene scores.
    #     # Will depend on the GWAS being in OT already.
    # }
    #### Set key and sort ####
    data.table::setkey(gene_hits, "SYMBOL")
    #### Save ####
    out_path <- save_snps2genes(
        gene_hits = gene_hits,
        save_path = save_path, 
        nCores = nCores,
        verbose = verbose
    )
    #### Report time ####
    report_time(start)
    #### Return ####
    if (isTRUE(return_path)) {
        return(out_path)
    } else {
        return(gene_hits)
    }
}
