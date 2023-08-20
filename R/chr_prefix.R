chr_prefix <- function(dat) {
    all(startsWith(as.character(dat$CHR[seq_len(10)]), "chr"))
}
