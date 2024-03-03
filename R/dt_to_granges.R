#' Convert data.table to GRanges object
#'
#' Convert \link[data.table]{data.table}/\link[base]{data.frame} to a 
#' \link[GenomicRanges]{GRanges} object.
#' 
#' @param dat Data.
#' @param chrom_col Chromosome column name.
#' @param start_col Genomic start position column name.
#' @param end_col Genomic end position column name.
#' @param style GRanges style (e.g. "NCBI, "UCSC") 
#' set by \link[GenomeInfoDb]{seqlevelsStyle}.
#' @param verbose Print messages.
#' 
#' @export
dt_to_granges <- function(dat,
                          chrom_col = "CHR",
                          start_col = "POS",
                          end_col = start_col,
                          style = "NCBI",
                          verbose = TRUE) {
    requireNamespace("GenomicRanges")
    requireNamespace("GenomeInfoDb")
    
    
    if (is_granges(dat)) {
        messager("dat is already a GRanges object.", v = verbose)
        gr.snp <- dat
    } else {
        messager("Converting dat to GRanges object.", v = verbose)
        dat <- data.table::as.data.table(dat)
        if(chrom_col=="seqnames"){
            data.table::setnames(dat,"seqnames","SEQnames")
        } 
        gr.snp <- GenomicRanges::makeGRangesFromDataFrame(
            dat,
            seqnames.field = chrom_col,
            start.field = start_col,
            end.field = end_col,
            keep.extra.columns = TRUE
        )
    }
    suppressMessages(suppressWarnings(
        GenomeInfoDb::seqlevelsStyle(gr.snp) <- style
    ))
    return(gr.snp)
}
