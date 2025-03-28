% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assignMutationsSignatureProbability.R
\name{assignSignatureProbabilityToMutations}
\alias{assignSignatureProbabilityToMutations}
\title{Assign Signature Probability to Mutations}
\usage{
assignSignatureProbabilityToMutations(
  sampleMutations,
  sampleSigsExposures,
  signatures,
  catalogue = NULL,
  enableUnassigned = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{sampleMutations}{list of annotated mutations obtained usually after building catalogues, for example using the functions vcfToSNVcatalogues or bedpeToRearrCatalogues}

\item{sampleSigsExposures}{matrix of exposures with only one row (one sample) and the exposures for each signature as columns}

\item{signatures}{mutational signatures matrix with the signatures used during signature fitting}

\item{catalogue}{sample catalogue, used only if enableUnassigned is true}

\item{enableUnassigned}{if true, probability can be assigned to an unassigned category, rather than a signature. This can lead to structural variants that are classified as unassigned. In practice, the positive part of the difference between catalogue and reconstruction is used as the unassigned signature and the unassigned exposures are also used.}
}
\value{
matrix of mutations with an additional column containing the probabilities in a text format that can be expanded into a matrix using expandColumnToMatrix
}
\description{
Given an annotated set of mutations and an exposure matrix from signature fit, this function estimates the probability that each mutation generates from each signature.
}
