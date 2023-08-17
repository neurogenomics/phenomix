chr_prefix <- function(dat) {
    all(startsWith(as.character(dat$CHR[seq(1, 10)]), "chr"))
}
