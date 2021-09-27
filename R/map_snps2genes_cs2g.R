#' Map SNPs to genes: cS2G
#'
#' Map SNPs (both coding and non-coding) onto genes
#' using the combined SNP-to-gene (cS2G) strategy..
#'
#' @section References:
#'
#' \href{https://doi.org/10.1101/2021.08.02.21261488}{
#' Gazal et al. 2021 preprint}
#'
#' \href{https://alkesgroup.broadinstitute.org/cS2G}{
#' Broad Institute server}
#'
#' @param sumstats_file GWAS summary statistics munged by
#' \link[MungeSumstats]{format_sumstats}.
#' Can be a path to the saved file or \link[data.table]{data.table}.
#' @param nCores Number of cores to parallelise across.
#' @param method Method to use for mapping SNPs to genes.
#' @param verbose Print messages.
#'
#' @return \code{gene_hits} \link[data.table]{data.table}
#'
#' @export
#' @importFrom data.table setkey fread
map_snps2genes_cs2g <- function(sumstats,
                                ref_genome = "GRCh37",
                                dataset = "Nasser2020",
                                abc = NULL,
                                nCores = 1,
                                verbose = TRUE) {

    #### Check build ####
    # if(!toupper(ref_genome) %in% c("HG19","GRCH37")){
    #     stop("sumstats must first be aligned to GRCh37, ",
    #          "as the ABC model results are aligned to GRCh37.")
    # }
    #### Import ABC predictions and merge with overlapping sumstats ####
    merged_hits <- merge_abc(
        sumstats = sumstats,
        dataset = dataset,
        abc = abc,
        nCores = nCores,
        verbose = verbose
    )
    #### Aggregate ABC ####
    gene_hits <- aggregate_sumstats(
        merged_hits = merged_hits,
        agg_var = "TargetGene",
        verbose = verbose
    )
    #### Get gene length ####
    gene_hits <- get_gene_length(
        gene_hits = gene_hits,
        ref_genome = ref_genome,
        gene_var = "TargetGene",
        use_symbols = TRUE,
        verbose = verbose
    )
    #### Compute Z-score ####
    zscore(
        dat = gene_hits,
        column = "BETA_mean",
        verbose = verbose
    )
    #### Rename TargetGene ####
    data.table::setnames(gene_hits, "TargetGene", "GENE")
    #### Return ####
    return(gene_hits)
}
