
#' HRD-LOH index from ASCAT
#' 
#' Compute the HRD-LOH index starting from a data fram (possibly loaded from an ASCAT output file)
#' with the following columns: 'seg_no', 'Chromosome', 'chromStart', 'chromEnd', 'total.copy.number.inNormal', 
#' 'minor.copy.number.inNormal', 'total.copy.number.inTumour', 'minor.copy.number.inTumour'
#' 
#' @param ascat.data data frame with the following columns: 'seg_no', 'Chromosome', 'chromStart', 'chromEnd', 'total.copy.number.inNormal', 'minor.copy.number.inNormal', 'total.copy.number.inTumour', 'minor.copy.number.inTumour'
#' @param SAMPLE.ID sample name
#' @return returns the HRD-LOH index
#' @references Abkevich, V., Timms, K. M., Hennessy, B. T., Potter, J., Carey, M. S., Meyer, L. a., ... Lanchbury, J. S. (2012). Patterns of genomic loss of heterozygosity predict homologous recombination repair defects in epithelial ovarian cancer. British Journal of Cancer, 107(10), 1776–82. https://doi.org/10.1038/bjc.2012.451
#' @export
#' @examples
#' ascat.data <- read.table("ascat.scv",sep=",",header=TRUE)
#' hrd_index <- ascatToHRDLOH(ascat.df,"test_sample")
ascatToHRDLOH <- function(ascat.data,SAMPLE.ID){

  #----------------------
  # Dominik's code BEGIN
  #----------------------
  # load ASCAT NGS profile
  # ascat.data <- read.table(ASCAT.PATH,header=FALSE, sep=',')
  # names(ascat.data ) <- c('seg_no', 'Chromosome', 'chromStart', 'chromEnd', 'total.copy.number.inNormal', 'minor.copy.number.inNormal', 'total.copy.number.inTumour', 'minor.copy.number.inTumour')
  ascat.data2 <- data.frame(SampleID=rep(SAMPLE.ID,nrow(ascat.data )),
                            Chromosome=ascat.data$Chromosome,
                            Start=ascat.data$chromStart,
                            End=ascat.data$chromEnd,
                            nProbes=NA,
                            totalCN=ascat.data$total.copy.number.inTumour,
                            nA=ascat.data$total.copy.number.inTumour - ascat.data$minor.copy.number.inTumour,
                            nB=ascat.data$minor.copy.number.inTumour,
                            Ploidy=rep(NA ,nrow(ascat.data )), #ploidy
                            AberrantCellFraction=rep(NA ,nrow(ascat.data )) #cellularity
  )
  
  ### reorder - HEADERS MUST BE IN CORRECT ORDER #####
  ll<-match(c("SampleID","Chromosome","Start","End","nProbes","totalCN","nA","nB","Ploidy" ,"AberrantCellFraction"),colnames(ascat.data2))
  ascat.data2 <-ascat.data2 [,ll]
  rm(ll)
  #----------------------
  # Dominik's code END
  #----------------------
  
  ascat.data2[,"Chromosome"] <- as.character(ascat.data2[,"Chromosome"])
  ascat.data2[ascat.data2[,"Chromosome"]=="X","Chromosome"] <- "23"
  ascat.data2[ascat.data2[,"Chromosome"]=="Y","Chromosome"] <- "24"
  ascat.data2[,"Chromosome"] <- as.numeric(ascat.data2[,"Chromosome"])
  
  HRD_LOH <- calc.hrd(ascat.data2, nA=7,check.names=FALSE, return.loc=FALSE)
  
  return(HRD_LOH)
}
