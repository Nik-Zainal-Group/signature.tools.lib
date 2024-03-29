% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SignatureFitLib.R
\name{exposureDistributionBarplot}
\alias{exposureDistributionBarplot}
\title{Distribution of Signatures in Samples}
\usage{
exposureDistributionBarplot(fileout = NULL, catalogue, exposures)
}
\arguments{
\item{fileout}{if specified, generate a plot, otherwise no plot is generated, use extension png or jpg}

\item{catalogue}{original catalogue, channels as rows and samples as columns}

\item{exposures}{exposures/activities of signatures in each sample. Signatures as rows, samples as columns}
}
\value{
list with to objects: the matrix of the distribution of the signatures in the samples and the hierarchical clustering object
}
\description{
Given a catalogue of samples and an exposures table, compute the relative amount of each signature
in each sample and the unassigned mutations. Also cluster the samples with hierarchical clustering
with average linkage and order the samples according to the clustering. Optionally, plot to file.
}
\examples{
res <- SignatureFit_withBootstrap(cat = catalogue,
                  signature_data_matrix = cosmic30,
                  nboot = 5,
                  threshold_percent = 0.1,
                  threshold_p.value = 0.1)
distribution_object <- exposureDistributionBarplot(catalogue=catalogue,
                  exposures=res$E_median_filtered)
}
\keyword{samples}
\keyword{unexplained}
