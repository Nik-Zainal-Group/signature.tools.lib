% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signatureFitMultiStepLib.R
\name{saveFitToFile}
\alias{saveFitToFile}
\title{Save fit object to file}
\usage{
saveFitToFile(fitObj, filename, verbose = T)
}
\arguments{
\item{fitObj}{object obtained from the Fit of FitMS function}

\item{filename}{file name where to save the fit object}
}
\description{
This function saves a Fit or FitMS object to an R data file using a standard name
}
\examples{
fitObj <- FitMS(catalogues,organ="Breast")
saveFitToFile(fitObj,"fit.rData")
}
