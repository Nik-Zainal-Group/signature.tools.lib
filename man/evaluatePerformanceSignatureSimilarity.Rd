% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performanceEvaluationLib.R
\name{evaluatePerformanceSignatureSimilarity}
\alias{evaluatePerformanceSignatureSimilarity}
\title{evaluate performance signature similarity}
\usage{
evaluatePerformanceSignatureSimilarity(
  true_signatures,
  estimated_signatures,
  true_exposures = NULL,
  outfile = NULL
)
}
\arguments{
\item{true_signatures}{matrix of true signatures, for example obtained using simulations, with mutational signatures names as colnames and channels as rownames}

\item{estimated_signatures}{matrix of estimated signatures, for example obtained using a signature extraction on a simulated dataset, with mutational signatures names as colnames and channels as rownames.}

\item{true_exposures}{matrix of true exposures, for example obtained using simulations, with mutational signatures names as colnames and sample names as rownames. This is optional, and it useful to show which of the true_signature are actually present in the dataset.}

\item{outfile}{file name for plotting, please use .pdf file name extension. Can be omitted.}
}
\value{
cosine similarity and signatures match between true and estimated signature names
}
\description{
Compare a signatures matrix of true signatures with a matrix of estimated signatures, and
compute a match between the signatures and their cosine similarity. The match implies
an optimal assignment that maximises the cosine similarity while allowing each true_signature
to match at most one estimated signature and viceversa.
If the number of signatures differ, the least similar signatures are removed until the number
of signatures is the same
Results are plotted to a file if outfile is specified.
}
