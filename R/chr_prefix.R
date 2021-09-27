chr_prefix <- function(sumstats) {
    all(startsWith(as.character(sumstats$CHR[seq(1, 10)]), "chr"))
}
