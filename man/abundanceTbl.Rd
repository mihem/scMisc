% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seurat_helper.R
\name{abundanceTbl}
\alias{abundanceTbl}
\title{wrapper function to calculate absolute and relative abundance}
\usage{
abundanceTbl(object, row_var, col_var, dir_output = ".")
}
\arguments{
\item{object}{Seurat object}

\item{row_var}{variable in meta data that will represent the rows}

\item{col_var}{variable in meta data that will represent the columns}

\item{dir_output}{target directory to save the results (default: .)}
}
\description{
wrapper function to calculate absolute and relative abundance
}
\examples{
library(Seurat)
abundanceTbl(pbmc_small, row_var = "groups", col_var = "letter.idents")
unlink("abundance_tbl_pbmc_small_letter.idents.xlsx")
}
