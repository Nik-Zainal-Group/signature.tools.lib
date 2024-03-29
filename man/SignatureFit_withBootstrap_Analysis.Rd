% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SignatureFitLib.R
\name{SignatureFit_withBootstrap_Analysis}
\alias{SignatureFit_withBootstrap_Analysis}
\title{Mutational Signatures Fit with Bootstrap Analysis}
\usage{
SignatureFit_withBootstrap_Analysis(
  outdir,
  cat,
  signature_data_matrix,
  nboot = 100,
  type_of_mutations = "subs",
  threshold_percent = 5,
  threshold_p.value = 0.05,
  method = "KLD",
  bf_method = "CosSim",
  alpha = -1,
  doRound = FALSE,
  nparallel = 1,
  n_sa_iter = 500
)
}
\arguments{
\item{outdir}{output directory for the analysis, remember to add '/' at the end}

\item{cat}{catalogue matrix, patients as columns, channels as rows}

\item{signature_data_matrix}{signatures, signatures as columns, channels as rows}

\item{nboot}{number of bootstraps to use, more bootstraps more accurate results}

\item{type_of_mutations}{either "subs", "rearr" or "generic"}

\item{threshold_percent}{threshold in percentage of total mutations in a sample, only exposures larger than threshold are considered}

\item{threshold_p.value}{p-value to determine whether an exposure is above the threshold_percent. In other words, this is the empirical probability that the exposure is lower than the threshold}

\item{method}{KLD or NNLS or SA}

\item{bf_method}{bleeding filter method, one of KLD or CosSim, only if bleeding filter is used (alpha>-1)}

\item{alpha}{set alpha to -1 to avoid Bleeding Filter}

\item{doRound}{round the exposures to the closest integer}

\item{nparallel}{to use parallel specify >1}

\item{n_sa_iter}{set max Simulated Annealing iterations if method==SA}

\item{verbose}{use FALSE to suppress messages}
}
\value{
returns the activities/exposures of the signatures in the given sample and other information, such as p-values and exposures of individual bootstrap runs.
}
\description{
This function is a wrapper for the function SignatureFit_withBootstrap_Analysis, which
produces several plots for each sample in the catalogues cat.
Fit a given set of mutational signatures into mutational catalogues to extimate the
activty/exposure of each of the given signatures in the catalogues. Implementation
of method similar to Huang 2017, Detecting presence of mutational signatures with
confidence, which uses a bootstrap apporach to calculate the empirical probability
of an exposure to be larger or equal to a given threshold (i.e. 5% of mutations of a sample).
This probability can be used to decide which exposures to remove from the initial fit,
thus increasing the sparsity of the exposures.
Note that SignatureFit_withBootstrap_Analysis will save the results of SignatureFit_withBootstrap
in the outdir directory using the R save() function. If SignatureFit_withBootstrap_Analysis
is rerun with the same setting, the saved file will be loaded to avoid rerunning the Signature Fit
and figures will be replotted.
}
\examples{
res <- SignatureFit_withBootstrap_Analysis(catalogues,signature_data_matrix)
}
\keyword{fit}
\keyword{mutational}
\keyword{signatures}
