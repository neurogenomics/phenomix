#' Map SNPs to genes: txdb
#' 
#' Map both coding and non-coding SNPs to genes using 
#' #' \link[TxDb.Hsapiens.UCSC.hg19.knownGene]{TxDb.Hsapiens.UCSC.hg19.knownGene}
#' or
#' \link[TxDb.Hsapiens.UCSC.hg38.knownGene]{TxDb.Hsapiens.UCSC.hg38.knownGene}
#' as a reference.
#' 
#' Permits selection of different region types with 
#' \code{queries} argument.
#' 
#' @param sumstats GWAS summary statistics munged by 
#' \link[MungeSumstats]{format_sumstats}. 
#' @param promoter_upstream How many basepairs upstream of known promoters
#'  to search for genes.
#' @param promoter_downstream How many basepairs dowmstream of known promoters
#' to search for genes.
#' @param queries Genomic region types to search for overlap with \code{gr}.
#' See \link[VariantAnnotation]{AllVariants} for details. 
#' @param all_variants Search all genomic regions regardless of type.
#' This can take longer and returns less targeted results.
#' @param verbose Print messages. 
#' @inheritParams MungeSumstats::format_sumstats
#' 
#' @return \code{hits} \link[data.table]{data.table}
#' 
#' @keywords internal
#' @importFrom VariantAnnotation locateVariants  
#' @importFrom VariantAnnotation PromoterVariants CodingVariants AllVariants IntronVariants ThreeUTRVariants SpliceSiteVariants
#' @importFrom GenomicRanges GRangesList
#' @importFrom parallel mclapply
map_snps_txdb <- function(sumstats,
                          promoter_upstream=35000,
                          promoter_downstream=10000,
                          queries=c("promoter","coding","intron",
                                    "threeUTR","splicesite"),
                          all_variants=FALSE,
                          ref_genome="GRCh37",
                          return_merged=TRUE,
                          nCores=1,
                          verbose=TRUE){  
    
    #### Convert to GRanges #### 
    gr <- MungeSumstats:::to_GRanges(sumstats)
    GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
    #### Get txdb ####
    txdb <- select_txdb_build(ref_genome = ref_genome) 
    ##### Gather region queries #####
    promoter_region <- VariantAnnotation::PromoterVariants(upstream=promoter_upstream, 
                                                           downstream=promoter_downstream)
    if(all_variants){
        # regions <- VariantAnnotation::AllVariants(promoter = promoter_region)
        regions <- list(allvariants=VariantAnnotation::AllVariants())
    }else {
        regions <- list(promoter=promoter_region, 
                        coding=VariantAnnotation::CodingVariants(),
                        intron=VariantAnnotation::IntronVariants(),
                        threeUTR=VariantAnnotation::ThreeUTRVariants(),
                        splicesite=VariantAnnotation::SpliceSiteVariants() 
        )
    } 
    messager(paste0("Using ",length(regions)," genomic region types:\n",
             paste0(" - ",names(regions),collapse = "\n")),
             v=verbose)
    regions <- regions[queries[queries %in% names(regions)]]
    #### Iterate over each region class ####
    hits <- parallel::mclapply(names(regions), function(x){
        message_parallel("Querying regions: ",x)
        res <- VariantAnnotation::locateVariants(query = gr, 
                                                 subject = txdb,
                                                 region = regions[[x]])
        res$region <- x
        return(res)
    }, mc.cores = nCores) %>% `names<-`(names(regions)) %>%
        GenomicRanges::GRangesList() %>%
        unlist() 
    #### Add RSID back in so we can merge later ####
    hits$SNP <- gr$SNP[hits$QUERYID] 
    
    #### Merge ####
    if(return_merged){
        merged_hits <- merge(x = sumstats, 
                             y = data.table::data.table(data.frame(hits)),
                             by="SNP", all=TRUE) 
        return(merged_hits)
    } else {
        return(hits)
    }  
}