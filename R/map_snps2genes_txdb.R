#' Map SNPs to genes: txdb
#'
#' Map SNPs (both coding and non-coding) onto genes
#' using
#' \link[TxDb.Hsapiens.UCSC.hg19.knownGene]{TxDb.Hsapiens.UCSC.hg19.knownGene}
#' or
#' \link[TxDb.Hsapiens.UCSC.hg38.knownGene]{TxDb.Hsapiens.UCSC.hg38.knownGene}
#' as a reference.
#'
#' @param dat GWAS summary statistics munged by
#' \link[MungeSumstats]{format_sumstats}.
#' Can be a path to the saved file or \link[data.table]{data.table}.
#' @param nCores Number of cores to parellelise across.
#' @param method Method to use for mapping SNPs to genes.
#' @param verbose Print messages.
#' @inheritParams map_snps_txdb
#' @inheritParams aggregate_sumstats 
#' @inheritParams translate_geneids_txdb
#'
#' @return \code{gene_hits} \link[data.table]{data.table}
#'
#' @keywords internal
#' @importFrom data.table setkey fread
map_snps2genes_txdb <- function(dat,
                                promoter_upstream = 35000,
                                promoter_downstream = 10000,
                                queries = c(
                                    "promoter", "coding", "intron",
                                    "threeUTR", "splicesite"
                                ),
                                drop_na = TRUE,
                                agg_var = "SYMBOL",
                                zscore_col = "BETA_mean",
                                nCores = 1,
                                verbose = TRUE) {
    # devoptera::args2vars(map_snps2genes_txdb)
    
    start <- Sys.time()
    #### Map SNPs to genes ####
    merged_hits <- map_snps_txdb(
        dat = dat,
        promoter_upstream = promoter_upstream,
        promoter_downstream = promoter_downstream,
        queries = queries,
        nCores = nCores
    )
    #### Translate gene IDs back to gene symbols ####
    merged_hits <- translate_geneids_txdb(
        gene_hits = merged_hits,
        drop_na = drop_na,
        verbose = verbose
    )
    #### Aggregate dat by gene #####
    gene_hits <- aggregate_sumstats(
        merged_hits = merged_hits,
        agg_var = agg_var,
        verbose = verbose
    ) 
    #### Get gene length ####
    gene_hits <- get_gene_length(
        gene_hits = gene_hits,
        drop_na = drop_na,
        verbose = verbose
    )
    #### Compute Z-score ####
    zscore(
        dat = gene_hits,
        column = zscore_col,
        verbose = verbose
    )
    #### Report time ####
    messager(utils::capture.output(difftime(Sys.time(), start)))
    ### Currently missing (relative to MAGMA):
    # P (tho we do have P_mean)
    # NPARAM
    #### Return ####
    return(gene_hits)
}
