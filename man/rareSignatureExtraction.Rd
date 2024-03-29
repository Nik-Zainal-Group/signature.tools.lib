% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rareSignaturesExtractionLib.R
\name{rareSignatureExtraction}
\alias{rareSignatureExtraction}
\title{Rare Signatures Extraction}
\usage{
rareSignatureExtraction(
  outfileRoot,
  catalogues,
  commonSignatures,
  commonExposures,
  residuals,
  unexpl_samples,
  clusters,
  useclusters,
  commonSigsToIgnore = NULL,
  checkForMissed = FALSE,
  maxiter = 1000,
  nsigs = 1
)
}
\arguments{
\item{outfileRoot}{if specified, generate a plot, otherwise no plot is generated}

\item{catalogues}{original catalogues, channels as rows and samples as columns}

\item{commonSignatures}{common mutational signatures used for fitting, channels as rows, signatures as columns}

\item{commonExposures}{exposures obtained from fitting catalogues using commonSignatures, signatures names as row names and sample names as column names. Typically this should be obtained from the unexplainedSamples function.}

\item{residuals}{residuals from catalogues fitted using common signatures. Typically this should be obtained from the unexplainedSamples function.}

\item{unexpl_samples}{names of samples that were found to be unexplained by the unexplainedSamples function}

\item{clusters}{vector of integers indicating the cluster that each unexplained sample belongs to. Typically this is obtained clustering residuals of unexplained samples using the catalogueClustering function.}

\item{useclusters}{list object indicating which cluster or group of clusters should be used to extract each rare signature. This is a list object and the length of the list indicates how many rare signatures will be extracted. Each entry of the list is a vector indicating one or more cluster numbers that will be combined into one for a rare signature extraction. For example, list(c(1,4),c(3),c(2)) indicates that three rare signatures will be extracted, one using samples from clusters 1 and 4, one using samples in cluster 3, and one using samples in cluster 2.}

\item{commonSigsToIgnore}{this parameter is used to reduce the importance of specific common signatures, by setting their exposures to 0 in the commonExposures table, thus changing the starting point of the optimisation. This is very useful when the rare signature to be extracted bears similarities with common signatures, which in turn interfere with the extraction. Leave NULL if not used, or use NA for extraction where this should not be used. For example, list(NA,c("commonSig1","commonSig3"),NA) assumes that three rare signatures will be extracted (determined by the useclusters parameter) and that the second extraction should ignore commonSig1 and commonSig3, which should also be row names of commonExposures.}

\item{checkForMissed}{if set to TRUE, there will be an additional check where each extracted rare signature will be compared to all the provided residuals to identify additional samples that may have the rare signature (min cosine similarity signature vs residual 0.95)}

\item{maxiter}{maximum number of iterations for the Lee and Seung multiplicative algorithm. Can be a single value or a vector if different extractions should have different maxiter values. It should be noted that performing an extraction on a single sample has potentially infinite solutions. For this reason it is useful in this case to run the Lee and Seung multiplicative algorithm starting from the commonExposures solution and run a limited number of iterations to find a neighbouring solution}

\item{nsigs}{How many signatures should be extracted. Typically this should be left to 1, indicating that for each group of unexplained samples (indicated by each entry in useclusters), only one rare signature should be extracted. If multiple rare signatures are thought to be present in one of the groups, then nsigs can be a vector (same length as useclusters) and a different values can be used, for example c(1,2,1) will extract 2 rare signatures from the samples in the second group.}
}
\value{
list of rare signatures and list of sample names for each signature
}
\description{
Extract rare signatures from sample catalogues. This function is the last step of
the extraction of rare signatures from sample catalogues. Before running this
function one should use the unexplainedSamples function to identify which
samples are likely to contain a rare signature. The unexplainedSamples function
produces a residual for each sample, indicating the part of the catalogue that is
not explained by the common signatures used. Then, one should cluster the residuals
of the samples that the unexplainedSamples function considered significant, using the
cataloguesClustering function. Finally, one should be ready to run this function.
See Examples below. This function will use the NNLM package to extract rare functions
while at the same time fitting the common signatures. An extraction is performed
for each of the clusters of residuals obtained using cataloguesClustering.
The extraction will use a variant of the Lee and Seung multiplicative algorithm,
optimising the KLD. The starting point for the extraction will be the commonExposures
provided. Specific signature exposures can be ignored with the parameter commonSigsToIgnore.
}
\examples{
resUnexpl <- unexplainedSamples(catalogues=catalogues,
                                sigs=signatures)
significant_residuals <- resUnexpl$all_residuals[,resUnexpl$which_significant]
clusterResiduals <- cataloguesClustering(significant_residuals,
                                         nclusters = 1:5)
resObj <- rareSignatureExtraction(catalogues=catalogues,
                                  commonSignatures=signatures,
                                  commonExposures=resUnexpl$exposures,
                                  residuals=resUnexpl$all_residuals,
                                  unexpl_samples=resUnexpl$unexplSamples,
                                  clusters=clusterResiduals$clusters_table[,"3"],
                                  useclusters=list=(c(1),c(3)),
                                  maxiter=c(100,1000))
}
\keyword{extract}
\keyword{rare}
\keyword{signatures}
