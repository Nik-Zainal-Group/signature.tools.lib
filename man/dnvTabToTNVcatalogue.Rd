% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TrinucleotideVariantsLib.R
\name{dnvTabToTNVcatalogue}
\alias{dnvTabToTNVcatalogue}
\title{Build a Triinucleotide Variants Catalogue from DNVs}
\usage{
dnvTabToTNVcatalogue(dnvtab)
}
\arguments{
\item{dnvtab}{requires columns Sample, Chrom, Pos, Ref, Alt, with Ref and Alt of length 2}
}
\value{
list of TNVs and TNV catalogue
}
\description{
This function takes as input a list of double nucleotide variants (DNVs),
and computes a list of trinucleotide variants (TNVs) finding which DNVs
are next to each other. It then returns the annotated TNVs and the TNV catalogues.
The trinucleotide mutations considered are only those where all three bases change,
i.e. TTT>AAA but not TTT>ATA.
}
