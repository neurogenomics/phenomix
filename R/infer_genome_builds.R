#' Infer genome builds
#' 
#' Infers the genome build of  summary statistics files (GRCh37 or GRCh38) 
#' from the data. Uses SNP (RSID) & CHR & BP to get genome build.
#'
#' @param nCores Number of threads to use for parallel processes.  
#' @param verbose Print messages.
#' @inheritParams MungeSumstats::get_genome_build
#' 
#' @return ref_genome the genome build of the data
#'
#' @examples
#' #Pass path to Educational Attainment Okbay sumstat file to a temp directory
#'
#' eduAttainOkbayPth <- system.file("extdata","eduAttainOkbay.txt", 
#'                                  package="MungeSumstats")
#' 
#' ## Call uses reference genome as default with more than 2GB of memory,
#' ## which is more than what 32-bit Windows can handle so remove certain checks
#' is_32bit_windows <- 
#' .Platform$OS.type == "windows" && .Platform$r_arch == "i386"
#' if (!is_32bit_windows) {
#' ref_genome <- infer_genome_builds(sumstats=eduAttainOkbayPth)                                  
#'}
#' @export
#' @importFrom MungeSumstats get_genome_build
#' @importFrom parallel mclapply
#' @importFrom utils capture.output
#' @importFrom dplyr %>%
infer_genome_builds <- function(sumstats,
                                header_only=TRUE,
                                sampled_snps=10000,
                                nCores=1,
                                verbose=TRUE){ 
    start <- Sys.time()
    sumstats <- unique(sumstats)
    messager("Inferring genome build of",length(sumstats),
             "sumstats file(s).",
             v=verbose)
    #### Get basename ####
    munged_ids <- gsub(paste(MungeSumstats:::supported_suffixes(),collapse = "|"),
                       "",basename(sumstats))
    #### Infer builds ####
    builds <- parallel::mclapply(sumstats, 
                                 function(x,
                                          .sampled_snps=sampled_snps){
        message_parallel(x)
        MungeSumstats::get_genome_build(sumstats = x, 
                                        sampled_snps = .sampled_snps,
                                        nThread = 1)
    }, mc.cores = nCores) %>% `names<-`(munged_ids) 
    #### Report time ####
    messager(utils::capture.output(difftime(Sys.time(), start)))
    #### Report build counts ####
    build_counts <- table(unlist(builds))
    for(i in seq(1,length(build_counts))){
        messager(names(build_counts)[i],":",build_counts[i],"file(s)",v=verbose)
    } 
    return(builds)
}
