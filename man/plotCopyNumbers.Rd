% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotCopyNumbers.R
\name{plotCopyNumbers}
\alias{plotCopyNumbers}
\title{plot Copy Numbers from ASCAT-like file}
\usage{
plotCopyNumbers(
  sv_df,
  sample_name,
  filename = NULL,
  plottitle = NULL,
  mar = NULL,
  minorCNcolour = "green",
  totalCNcolour = "red",
  outofrangeCNcolour = "purple",
  highlightRegions = NULL,
  highlightText = NULL,
  highlightColour = "brown",
  textscaling = 1,
  genome.v = "hg19"
)
}
\arguments{
\item{sv_df}{data frame with the following columns: 'seg_no', 'Chromosome', 'chromStart', 'chromEnd', 'total.copy.number.inNormal', 'minor.copy.number.inNormal', 'total.copy.number.inTumour', 'minor.copy.number.inTumour'}

\item{sample_name}{sample name}

\item{filename}{if a file name is specified, a pdf file will be generated}

\item{mar}{set custom margins, otherwise mar will be set to mar = c(4,4,3,2)}

\item{highlightRegions}{dataframe with columns 'Chromosome', 'Start', 'End', which specify the regions to highlight. This is optional.}

\item{highlightText}{if highlightRegions is specified, then highlightText will be plotted in the bottom right}

\item{highlightColour}{set the colour of the highlight regions and text. Default is brown}
}
\value{
returns the HRD-LOH index or regions
}
\description{
Plot the copy numbers across the Chromosomes. Optionally, plot also highlight regions at the bottom. For example, could be used
to highlight HRD-LOH regions.
}
\examples{
ascat.data <- read.table("ascat.scv",sep=",",header=TRUE)
HRDLOHregions <- ascatToHRDLOH(ascat.df,"test_sample",return.loc=TRUE)
HRDLOHscore <- nrow(HRDLOHregions)
plotCopyNumbers(sv_df = ascat.data,sample_name = "",highlightRegions = HRDLOHregions,
                highlightText = paste0("HRD-LOH\\\\n",HRDLOHscore))
}
