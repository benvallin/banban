% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_n0.R
\name{calculate_n0}
\alias{calculate_n0}
\title{Calculate mean n0 from LinRegPCR analysis}
\usage{
calculate_n0(path, linregpcr_files, type, include_all = FALSE, hkg = NULL)
}
\arguments{
\item{path}{a string specifying the path to the directory containing LinRegPCR output files.}

\item{linregpcr_files}{a character vector containing the names of LinRegPCR output files, with the extension included (.xls, .xlsx or .csv).}

\item{type}{a string specifying the extension of LinRegPCR output files (either .xls, .xlsx or .csv).}

\item{include_all}{a logical, default is FALSE. If TRUE, will also include the samples that did not reach the plateau phase or displayed deviating PCR efficiency.}

\item{hkg}{optional. A named list of character vectors specifying the sets of amplicons to be used for computing reference mean n0 values. For each vector of the list, reference mean n0 values are named after the vector and computed as the geometric mean of the vector elements.}
}
\value{
A tibble containing the mean n0 value of each amplicon (including newly created references) for each RT sample.
}
\description{
calculate_n0() computes the mean n0 values of qPCR replicates from xls/xlsx/csv files containing LinRegPCR output.
}
\examples{

}
