% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{fPlotCustom}
\alias{fPlotCustom}
\title{Seurat feature plot}
\usage{
fPlotCustom(
  object,
  markers,
  par,
  reduction,
  width = 16,
  height = ceiling(length(genes_found)/4) * 3,
  dir_output = "."
)
}
\arguments{
\item{object}{Seurat object}

\item{markers}{a data frame with a column called \code{cell_source} that represents the cell population and its source and a column \code{gene}}

\item{par}{a character string representing the cell_source to plot}

\item{reduction}{a character string specifying the dimension reduction}

\item{width}{width of output plot (default: 16)}

\item{height}{height of output plot (default: length of genes divided by four, ceiling, times three)}

\item{dir_output}{directory to save the output plot (default: ".")}
}
\value{
save feature plot
}
\description{
create and save a Seurat feature plot
}
\examples{
library(Seurat)
markers <- data.frame(cell_source = c("B", "B"), gene = c("MS4A1", "CD79A"))
fPlotCustom(
  object = pbmc_small,
  markers = markers,
  par = "B",
  reduction = "tsne"
)
unlink("fp_pbmc_small_B.png")
}
