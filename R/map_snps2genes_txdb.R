#' Map SNPs to genes: txdb
#' 
#' Map SNPs (both coding and non-coding) onto genes
#' using
#' \link[TxDb.Hsapiens.UCSC.hg19.knownGene]{TxDb.Hsapiens.UCSC.hg19.knownGene}
#' or
#' \link[TxDb.Hsapiens.UCSC.hg38.knownGene]{TxDb.Hsapiens.UCSC.hg38.knownGene}
#' as a reference.
#'  
#' @param sumstats_file GWAS summary statistics munged by 
#' \link[MungeSumstats]{format_sumstats}.
#' Can be a path to the saved file or \link[data.table]{data.table}.
#' @param nCores Number of cores to parellelise across. 
#' @param method Method to use for mapping SNPs to genes. 
#' @param verbose Print messages.
#' @inheritParams map_snps_txdb
#' 
#' @return \code{gene_hits} \link[data.table]{data.table}
#' 
#' @export
#' @importFrom data.table setkey fread
map_snps2genes_txdb <- function(sumstats, 
                                promoter_upstream=35000,
                                promoter_downstream=10000,
                                queries=c("promoter","coding","intron",
                                          "threeUTR","splicesite"),
                                nCores=1,
                                verbose=TRUE){ 
    
    start <- Sys.time() 
    #### Map SNPs to genes ####
    merged_hits <- map_snps_txdb(sumstats = sumstats,
                                 promoter_upstream = promoter_upstream, 
                                 promoter_downstream = promoter_downstream,
                                 queries = queries,
                                 nCores = nCores)   
    #### Aggregate sumstats by gene #####   
    gene_hits <- aggregate_sumstats(merged_hits = merged_hits,
                                    agg_var = "GENEID",
                                    verbose = verbose)
    #### Translate gene IDs back to gene symbols ####
    gene_hits <- translate_geneids_txdb(gene_hits = gene_hits, 
                                        verbose = verbose) 
    #### Get gene length ####
    gene_hits <- get_gene_length(gene_hits = gene_hits,
                                 verbose = verbose)
    #### Compute Z-score ####
    zscore(dat = gene_hits, 
           column = "BETA_mean",
           verbose = verbose)  
    #### Report time ####
    messager(utils::capture.output(difftime(Sys.time(), start)))
    ### Currently missing (relative to MAGMA):
    # P (tho we do have P_mean)
    # NPARAM   
    #### Return ####
    return(gene_hits)
}