% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_trait_cor.R
\name{plot_trait_cor}
\alias{plot_trait_cor}
\title{Plot top trait-trait correlations}
\usage{
plot_trait_cor(
  knn,
  top_n = 3,
  non_self = TRUE,
  group_var = "trait1",
  show_plot = TRUE,
  type = "heat"
)
}
\arguments{
\item{knn}{A melted similarity matrix or K-nearest neighbor graph.}

\item{top_n}{The number of top correlations to plot.}

\item{non_self}{Remove trait2 that are present in trait1.}

\item{group_var}{Column in \code{knn} data to group by 
when computing \code{top_n} traits.}

\item{show_plot}{Whether to print the plot.}

\item{type}{Plot type ("heat" or "bar").}
}
\value{
\code{ggplot} object.
}
\description{
Takes the output of \link[phenomix]{find_neighbors} as input to \code{knn}.
}
\examples{
obj <- get_HPO()[seq(100),]
knn <- find_neighbors(obj = obj,
                      var1_search = "cardio",
                      label_col = "HPO_label")
gg_cor <- plot_trait_cor(knn=knn,  top_n = 3)
}
