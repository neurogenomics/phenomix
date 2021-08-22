select_txdb_build <- function(ref_genome,
                              verbose=TRUE){
    
    if(toupper(ref_genome) %in% c("HG19","GRCH37")){
        messager("Selecting TxDb.Hsapiens.UCSC.hg19.knownGene.",v=verbose)
        txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene    
    }else if(toupper(ref_genome)=="GRCH38") {
        messager("Selecting TxDb.Hsapiens.UCSC.hg38.knownGene",v=verbose)
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    } else {
        stop("ref_genome must be 'GRCh37' or 'GRCh38'")
    }
    return(txdb)
}