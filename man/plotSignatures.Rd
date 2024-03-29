% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signatureFitMultiStepLib.R
\name{plotSignatures}
\alias{plotSignatures}
\title{Plot Signatures with automated detection of type of mutations}
\usage{
plotSignatures(
  signature_data_matrix,
  output_file = NULL,
  plot_sum = TRUE,
  overall_title = "",
  add_to_titles = NULL,
  mar = NULL,
  howManyInOnePage = 100,
  ncolumns = 3
)
}
\arguments{
\item{signature_data_matrix}{matrix of signatures, signatures as columns and channels as rows}

\item{output_file}{set output file, should end with ".jpg", "png" or ".pdf". If output_file==null, output will not be to a file, but will still run the plot functions. The option output_file==null can be used to add this plot to a larger output file.}

\item{plot_sum}{whether the sum of the channels should be plotted. If plotting signatures this should be FALSE, but if plotting sample catalogues, this can be set to TRUE to display the number of mutations in each sample.}

\item{overall_title}{set the overall title of the plot}

\item{mar}{set the option par(mar=mar)}

\item{howManyInOnePage}{how many signatures or catalogues should be plotted on one page. Multiple pages are plotted if more signatures/catalogues to plot have been requested}

\item{ncolumns}{how many columns should be used to arrange the signatures/catalogues to plot}
}
\description{
This function checks the channels of the input matrix, determines the type of mutations and plots using the
most appropriate signature plot.
}
