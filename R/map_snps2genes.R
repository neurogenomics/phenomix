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
#' @param sumstats_file GWAS summary statistics munged by
#' \link[MungeSumstats]{format_sumstats}.
#' Can be a path to the saved file or \link[data.table]{data.table}.
#' @param adjust_z Whether to adjust Z-statistic using
#'  the \code{model} and \code{formula}.
#' @param nCores Number of cores to parallelise across.
#' @param method Method to use for mapping SNPs to genes.
#' @param save_dir Where to save results.
#' @param verbose Print messages.
#' @inheritParams import_abc
#' @inheritParams adjust_zstat
#'
#' @return gene_hits \link[data.table]{data.table}.
#'
#' @export
#' @importFrom data.table setkey fread
map_snps2genes <- function(sumstats_file,
                           #### Methods args ####
                           method = c("abc", "txdb"),
                           dataset = "Nasser2020",
                           abc = NULL,
                           #### Adjust Z-stat args ####
                           adjust_z = TRUE,
                           drop_MHC = TRUE,
                           model = NULL,
                           log_vars = c("NSNPS", "NPARAM", "GENELEN"),
                           formula = ZSTAT ~ NSNPS + logNSNPS + NPARAM +
                               logNPARAM + GENELEN + logGENELEN,
                           #### Util args ####
                           nCores = 1,
                           save_dir = tempdir(),
                           return_path = FALSE,
                           verbose = TRUE) {

    # save_dir <- "/Volumes/bms20/projects/neurogenomics-lab/live/GWAS_sumstats/OpenGWAS"
    # metagwas_all <- read.csv(file.path(save_dir,"OpenGWAS_metadata.csv"), row.names = 1)
    # sumstats_file <- metagwas_all$path[1]
    # method="abc"; nCores <- 10; verbose=TRUE; dataset="Nasser2020";  ref_genome="GRCh37"; #abc=NULL;
    # model=NULL; drop_MHC=TRUE; log_vars=c("NSNPS","NPARAM","GENELEN");   formula=ZSTAT ~ NSNPS + logNSNPS + NPARAM + logNPARAM + GENELEN + logGENELEN

    start <- Sys.time()
    #### Check method ####
    method <- tolower(method[1])
    #### Prepare sumstats ####
    # Replace with MungeSumstats:::format_sumstats() once it's fixed.
    if (methods::is(sumstats_file, "character")) {
        messager("Importing sumstats_file.", v = verbose)
        sumstats <- data.table::fread(sumstats_file,
            nThread = nCores
        )
    } else {
        sumstats <- sumstats_file
    }
    #### Select method and map SNPs ####
    #### txdb ####
    if (method == "txdb") {
        #### Map ####
        gene_hits <- map_snps2genes_txdb(
            sumstats = sumstats,
            verbose = verbose
        )
        #### Adjust ZSTAT ####
        if (adjust_z) {
            gene_hits <- adjust_zstat(
                dat = gene_hits,
                model = model,
                log_vars = log_vars,
                drop_MHC = drop_MHC,
                formula = formula,
                verbose = verbose
            )
        }
    }
    #### ABC ####
    if (method == "abc") {
        #### Map ####
        gene_hits <- map_snps2genes_abc(
            sumstats = sumstats,
            abc = abc,
            dataset = dataset,
            nCores = nCores,
            verbose = verbose
        )
        #### Adjust ZSTAT ####
        if (adjust_z) {
            gene_hits <- adjust_zstat(
                dat = gene_hits,
                model = model,
                log_vars = log_vars,
                drop_MHC = drop_MHC,
                formula = formula,
                verbose = verbose
            )
        }
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
        save_dir = save_dir,
        sumstats_file = sumstats_file,
        method = method,
        nCores = nCores,
        verbose = verbose
    )
    #### Report time ####
    report_time(start)
    #### Return ####
    if (return_path) {
        return(out_path)
    } else {
        return(gene_hits)
    }
}
