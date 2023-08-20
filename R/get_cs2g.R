#' Import precomputed cS2G results
#' 
#' Precomputed results from the combined SNP-to-gene mapping model 
#' (\href{https://doi.org/10.1038/s41588-022-01087-y}{
#' Gazal et al., Nature Genetics, 2022}).\cr
#' 
#' \strong{gwas_catalog_cS2G}:\cr\cr
#' This compressed file contains cS2G score for 78,499 potentially causal
#'  SNP-gene-disease triplets from the NHGRI-EBI GWAS catalog
#'   (accessed file: gwas_catalog_v1.0-associations_e100_r2021-02-25.txt)\cr
#' PUBMEDID: PUBMED ID\cr
#' DISEASE.TRAIT: PHENOTYPE\cr
#' CHR: Chromosome\cr
#' POS: Position (hg19)\cr
#' SNP: SNP ID\cr
#' GENE: Gene ID\cr
#' cS2G: cS2G score\cr
#' INFO: S2G scores for the 10 main functional S2G
#'  (including the 7 cS2G constituent strategies)\cr
#' 
#' \strong{finemapping_cS2G_UKBB}:\cr\cr
#' Import 47 GWAS (from UK Biobank) mapped to genes using the
#' cS2G strategy.\cr\cr
#' This compressed file contains cS2G score for 138,716 potentially causal 
#' SNP-gene-disease triplets inferred from causal SNP-disease pairs with 
#' PIP>0.05 from functionally informed fine-mapping of 49
#'  UK Biobank diseases/traits.\cr
#' CHR: Chromosome\cr
#' BP: Position (hg19\cr
#' SNP: SNP ID\cr
#' DISEASE.TRAIT: PHENOTYPE\cr
#' PIP: Posterior inclusion probability\cr
#' GENE: Gene ID\cr
#' cS2G: cS2G score\cr
#' INFO: S2G scores for the 10 main functional S2G (including the 7 cS2G 
#' constituent strategies)\cr
#' This compressed file also contains a beta2 directory with PolyFun outputs
#'  for each traits (used in the omnigenicity analyses)\cr
#'  
#' @section References:
#' \href{https://doi.org/10.1101/2021.08.02.21261488}{
#' Gazal et al. 2021 preprint}
#' \href{https://zenodo.org/record/7754032}{Zenodo} 
#' @param dataset Dataset from the association 
#' \href{https://zenodo.org/record/7754032}{Zenodo} repository to use.
#' @param value_var Which column you want to fill the matrix with. 
#' Only used when \code{as_matrix==TRUE}.
#' @param min_value Minimum value of \code{value_var} to include in data.
#' Helpful when you're using a machine that can't handle large data.
#' @param save_dir Directory to cache dataset in.
#' @param as_granges Convert the data object t0 
#' \link[GenomicRanges]{GRanges} format.
#' @inheritParams phenomix_merge 
#' @inheritParams downloadR::downloader
#' @inheritParams data.table::dcast.data.table
#' @returns Named with two items: "data" and "obs" (observation metadata)
#'
#' @export
#' @import data.table
#' @import scKirby
#' @import echotabix
#' @importFrom tools R_user_dir
#' @importFrom downloadR downloader
#' @importFrom utils unzip 
#' @examples
#' meta <- get_cs2g_ukb()
get_cs2g <- function(dataset=c("finemapping_cS2G_UKBB",
                               "gwas_catalog_cS2G"),
    save_dir = tools::R_user_dir(package = "phenomix",
                                 which = "cache"),
    as_matrix = FALSE,
    as_granges = FALSE,
    value_var = c("cS2G","PIP"),
    min_value=0,
    formula = "gene ~ ID",
    nThread = 1,
    verbose = TRUE) {
    # devoptera::args2vars(get_cs2g, reassign = TRUE)
    
    requireNamespace("downloadR")
    value_var <- value_var[1] 
    dataset <- dataset[1] 
    #### Check args ####
    if(sum(c(as_matrix,as_granges))==2){
        messager(
            "WARNING: Only one of as_matrix and as_granges",
            "can be set to TRUE at once. Setting as_matrix=FALSE",v=verbose)
        as_matrix <- FALSE
    }
    #### Get URL ####
    host <- "https://zenodo.org/record/7754032/"
    URL <- if(dataset=="finemapping_cS2G_UKBB"){
        paste0(host,"finemapping_cS2G_UKBB.zip?download=1")
    } else if (dataset=="gwas_catalog_cS2G"){
        paste0(host,"gwas_catalog_cS2G.zip?download=1")
    }
    #### Download ####
    tmp <- downloadR::downloader(input_url = URL, 
                                 output_dir = save_dir, 
                                 nThread = nThread,
                                 verbose = verbose)
    tmp <- utils::unzip(tmp, 
                        exdir = gsub("\\.zip|\\?download=1","",tmp)) 
    dat <- data.table::fread(tmp[[1]])
    #### Parse ####
    if(dataset=="finemapping_cS2G_UKBB"){
        out <- get_cs2g_ukb(dat)
    } else if (dataset=="gwas_catalog_cS2G"){
        out <- get_cs2g_gwascatalog(dat)
    }
    dat <- out$data
    obs <- out$obs
    #### Check value_var exists ####
    if(!value_var %in% names(dat)){
        stopper(paste0("value_var=",shQuote(value_var)),
                "not present in data.")
    } 
    #### Filter ####
    if(min_value>0){ 
        messager("Filtering",paste0(shQuote(value_var),">=",min_value),
                 v=verbose)    
        dat <- dat[dat[[value_var]]>min_value,]
    }
    if (isTRUE(as_matrix)) {
        messager("Converting to matrix format.", v = verbose) 
        X <- data.table::dcast.data.table(
            data = dat,
            formula = formula,
            value.var = value_var,
            fill = 0,
            fun.aggregate = "mean"
        ) |> scKirby::to_sparse(verbose = verbose) 
        #### Ensure rownames match ####
        obs <- obs[colnames(X), ] 
        return(list(
            data = X,
            obs = obs
        ))
    } else {
        if(isTRUE(as_granges)){
            query_granges <- echotabix::construct_query(
                query_dat = dat, 
                query_start_col = "BP")
            return(list(data=query_granges,
                        obs=obs))
        } 
        return(list(
            data = dat,
            obs = obs
        ))
    }
}
