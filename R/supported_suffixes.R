supported_suffixes <- function (tabular = TRUE, 
                                tabular_compressed = TRUE, 
                                vcf = TRUE, 
                                vcf_compressed = TRUE) {
    supported <- c()
    suffixes <- c(".tsv", ".txt", ".csv")
    suffixes.gz <- c(paste0(suffixes, ".gz"), paste0(suffixes, 
                                                     ".bgz"))
    suffixes.vcf <- c(".vcf")
    suffixes.vcf.gz <- c(paste0(suffixes.vcf, ".gz"), paste0(suffixes.vcf, 
                                                             ".bgz"))
    if (tabular) 
        supported <- c(supported, suffixes)
    if (tabular_compressed) 
        supported <- c(supported, suffixes.gz)
    if (vcf) 
        supported <- c(supported, suffixes.vcf)
    if (vcf_compressed) 
        supported <- c(supported, suffixes.vcf.gz)
    return(supported)
}
