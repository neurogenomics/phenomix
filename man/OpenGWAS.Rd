% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{OpenGWAS}
\alias{OpenGWAS}
\title{OpenGWAS metadata}
\format{
An object of class \code{data.table} (inherits from \code{data.frame}) with 42349 rows and 31 columns.
}
\usage{
OpenGWAS
}
\description{
Static copy of metadata for all GWAS/QTL datasets on
 \href{https://gwas.mrcieu.ac.uk/}{OpenGWAS}.
}
\details{
\code{
OpenGWAS <- MungeSumstats::find_sumstats() 
OpenGWAS$query = NULL
OpenGWAS <- opengwas_meta(meta=OpenGWAS)
OpenGWAS[,magma_done:=RCurl::url.exists(magma_out)]
usethis::use_data(OpenGWAS, overwrite = TRUE)
}
}
\keyword{datasets}
