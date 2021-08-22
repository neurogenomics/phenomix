#' Get MHC genes
#' 
#' Retrieve genes from the MHC 
#' (Major Histocompatibility Complex) region.
#' 
#' The MHC region tends to be ubiquitously associated with all GWAS,
#' so removing it allows us to focus on non-MHC signals of interest.
#' For further information related to the MHC region, see the sources below.
#' 
#' @source \href{https://www.nature.com/articles/s41467-019-11953-9}{DeGAs}
#'  
#' @param seqnames Genomic ranges to search for genes. 
#' @param force_new Query for all genes within \code{seqnames} using 
#' \link[TxDb.Hsapiens.UCSC.hg19.knownGene]{TxDb.Hsapiens.UCSC.hg19.knownGene}
#' or
#' \link[TxDb.Hsapiens.UCSC.hg38.knownGene]{TxDb.Hsapiens.UCSC.hg38.knownGene}
#' as a reference.
#' If \code{FALSE} \code{(DEFAULT)}, will instead 
#' return a stored \link[data.table]{data.table} 
#' aligned to the "GRCh37" reference genome.
#' @param verbose Print messages.
#' @inheritParams map_snps_txdb
#'  
#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom data.table data.table
#' 
#' @examples 
#'  MHC_genes <- get_mhc_genes()
get_mhc_genes <- function(seqnames="chr6:25000000-34000000",
                          ref_genome="GRCh37",
                          force_new=FALSE,
                          verbose=TRUE){
    
    if(force_new){
        messager("Retrieving MHC region genes with new query.",
                 v=verbose)
        gr_mhc <- GenomicRanges::GRanges(seqnames = seqnames)
        MHC_genes <- map_snps_txdb(gr = gr_mhc, 
                                 ref_genome = ref_genome,
                                 all_variants = TRUE) 
        MHC_genes <- data.table::data.table(data.frame(MHC_genes))
        MHC_genes <- translate_geneids_txdb(gene_hits = MHC_genes)
    } else {
        MHC_genes <- phenomix::MHC_genes
    } 
    return(MHC_genes)
}
