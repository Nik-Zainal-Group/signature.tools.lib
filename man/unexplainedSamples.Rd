% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SignatureFitLib.R
\name{unexplainedSamples}
\alias{unexplainedSamples}
\title{Estimate samples not fully explained by signature fit}
\usage{
unexplainedSamples(fileout = NULL, catalogue, sigs, exposures,
  pvalue_threshold = 0.01)
}
\arguments{
\item{fileout}{if specified, generate a plot, otherwise no plot is generated}

\item{catalogue}{original catalogue, channels as rows and samples as columns}

\item{sigs}{mutational signautures used for fitting, channels as rows, signatures as columns}

\item{exposures}{exposures/activities of signatures in each sample. Signatures as rows, samples as columns}

\item{pvalue_threshold}{threshold for statistical significance}
}
\value{
table of unexplained samples
}
\description{
Given a catalogue of samples, signatures and exposures, compute the sum of the absolute deviations (SAD)
between the original catalogue and the reconstructed samples (i.e. signatures x exposures) and
normalise this sum by the total number of mutations in the sample. Then, for each sample, compare its
normalised SAD to the normalised SAD of the other samples and check if it is significantly different. In practice,
a p-value is computed fitting a gaussian distribution to the other samples.
}
\examples{
res <- SignatureFit_withBootstrap(cat = catalogue,
                  signature_data_matrix = cosmic30,
                  nboot = 5,
                  threshold_percent = 0.1,
                  threshold_p.value = 0.1)
s_table <- unexplainedSamples(catalogue=catalogue,
                  sigs=cosmic30,
                  exposures=res$E_median_filtered)
}
\keyword{samples}
\keyword{unexplained}
