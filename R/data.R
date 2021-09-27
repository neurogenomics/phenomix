#' MHC genes
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
#' \code{
#' MHC_genes <- get_mhc_genes()
#' usethis::use_data(MHC_genes, overwrite = TRUE)
#' }
"MHC_genes"
