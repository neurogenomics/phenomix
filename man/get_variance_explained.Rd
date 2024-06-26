% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_variance_explained.R
\name{get_variance_explained}
\alias{get_variance_explained}
\title{Get variance explained}
\usage{
get_variance_explained(
  obj,
  reduction = names(obj@reductions)[1],
  layer = "scale.data",
  dims = NULL
)
}
\description{
Get the proportion of total variance explained
 by a given dimensionality reduction 
 (0-1 where 1 indicates 100% of variance explained).
}
\examples{
obj <- get_HPO()
out <- get_variance_explained(obj)
}
