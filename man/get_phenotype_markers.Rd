% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_trait_markers.R
\name{get_phenotype_markers}
\alias{get_phenotype_markers}
\title{Identify marker features for each phenotype}
\usage{
get_phenotype_markers(
  seurat,
  terms = NULL,
  assay1 = "MAGMA",
  assay1_slot = "scale.data",
  assay2 = "specificity",
  assay2_slot = "counts",
  n_features = 100
)
}
\arguments{
\item{seurat}{\code{Seurat} object.}

\item{terms}{A list of substrings to search for in column names.}

\item{assay1}{First assay name.}

\item{assay1_slot}{Slot to use from \code{assay1}.}

\item{assay2}{Second assay name.}

\item{assay2_slot}{Slot to use from \code{assay2}.}

\item{n_features}{Number of features to select per assay.}
}
\description{
Identify marker features for a set of phenotypes
by taking the intersection between
the top \code{n_features} in two assays.
}
