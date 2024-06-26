% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_cluster_colors.R
\name{add_cluster_colors}
\alias{add_cluster_colors}
\title{Add cluster colors}
\usage{
add_cluster_colors(
  obj,
  cluster_col = "seurat_clusters",
  color_col = "cluster_colors",
  preferred_palettes = "kovesi.cyclic_mygbm_30_95_c78",
  force_new = FALSE,
  return_dict = FALSE
)
}
\description{
Add or extract cluster colors from a \link{Seurat} object.
}
\examples{
obj <- get_HPO()
obj2 <- add_cluster_colors(obj)
}
