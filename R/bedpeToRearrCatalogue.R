# source('../lib/utils/rearrangement.clustering_bedpe.R')
# source('../lib/utils/prepare.rearr.catalogue_fromAnnotatedBedpe.R')
# source('../lib/classifyRearrangementsFromBedpe.R')


#' bedpeToRearrCatalogue
#'
#' This function converts a data frame BEDPE into a rearrangement catalogue, 
#' you should pass rearrangements of only one sample.
#' 
#' Columns present in the BEDPE should be:
#' 
#' "chrom1"
#' 
#' "start1"
#' 
#' "end1"
#' 
#' "chrom2"
#' 
#' "start2"
#' 
#' "end2"
#' 
#' "sample" sample name
#' 
#' "name" name of mate1
#' 
#' "partner" name of mate2
#' 
#' "score"
#' 
#' In addition, either the strands of the mates:
#' 
#' "strand1": + or -
#' 
#' "strand2": + or -
#' 
#' or the structural variant class
#' 
#' "svclass": translocation, inversion, deletion, tandem-duplication
#' 
#' According to Sanger notation you should have sv class equal to: (strand1/strand2)
#' 
#' inversion (+/-), if mates on the same chromosome
#' 
#' inversion (-/+), if mates on the same chromosome
#' 
#' deletion (+/+), if mates on the same chromosome
#' 
#' tandem-duplication (-/-), if mates on the same chromosome
#' 
#' translocation, if mates are on different chromosomes
#' 
#' @param sv_bedpe data frame BEDPE as described above
#' @return returns the rearrangement catalogue for the given sample
#' @keywords bedpe, rearrangement
#' @export
#' @examples
#' vcf_sv_file.bedpe <- "sample.bedpe"
#' sv_bedpe <- read.table(vcf_sv_file.bedpe,sep = "\t",header = TRUE,
#'                      stringsAsFactors = FALSE,check.names = FALSE)
#' #build a catalogue from the bedpe file
#' res.cat <- bedpeToRearrCatalogue(sv_bedpe)
bedpeToRearrCatalogue <- function(sv_bedpe){
  
  #Annotate the bedpe if necessary
  
  #check whether column is.clustered is present,
  #if not, compute it
  if (! "is.clustered" %in% colnames(sv_bedpe)){
    clustering.result <- rearrangement.clustering_bedpe(sv_bedpe,
                                                        plot.path = NA,
                                                        kmin=10,                                                  
                                                        kmin.samples=1,
                                                        gamma.sdev=25,
                                                        PEAK.FACTOR=10,
                                                        thresh.dist=NA)
    sv_bedpe <- clustering.result$sv_bedpe
  }
  #check whether column svclass is present,
  #if not, compute it
  if (! "svclass" %in% colnames(sv_bedpe)){
    if ("strand1" %in% colnames(sv_bedpe) & "strand2" %in% colnames(sv_bedpe)){
      sv_bedpe <- classifyRearrangementsFromBedpe(sv_bedpe)
    }else{
      message("cannot classify rearrangements: svclass column missing, and cannot compute it because strand1 and strand2 are missing.")
    }
  }
  
  #now compute the catalogue
  prepare.rearr.catalogue_fromAnnotatedBedpe(sv_bedpe)
}



