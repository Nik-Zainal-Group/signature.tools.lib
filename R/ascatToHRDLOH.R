#setwd("~/sandbox/HRDetectScripts/bin")

# source("../lib/hrdIndex/Functions.GenomicScars.r")
# 
# ASCAT.PATH <- NULL
# project_id <- NULL
# SAMPLE.ID <- NULL
# ploidy <- NULL
# cellularity <- NULL	
# work_dir_cn <- NULL
#   
# args = commandArgs(trailingOnly=TRUE)
# if (length(args)!=6) {
#   message("--------------------------------------")
#   message("I was expecting six inputs, but a number of inputs different from six was supplied. Entering debug mode!")
#   project_id <- "BLCA-US"
#   SAMPLE.ID <- "0c7aca3f-e006-4de3-afc2-20b4f727d4fd"
#   ploidy <- "2"
#   cellularity <- "0.6"	
#   work_dir_cn <- paste0("../results/pancan/",project_id,"/cn/")
#   ASCAT.PATH <- paste0(work_dir_cn,'0c7aca3f-e006-4de3-afc2-20b4f727d4fd.ascat_ngs.summary.csv')
#   message("ENTERING DEBUG MODE, using the following settings:")
#   message("ASCAT.PATH ",ASCAT.PATH)
#   message("project_id ",project_id)
#   message("SAMPLE.ID ",SAMPLE.ID)
#   message("ploidy ",ploidy)
#   message("cellularity ",cellularity)
#   message("work_dir_cn ",work_dir_cn)
#   message("--------------------------------------")
# } else {
#   ASCAT.PATH <- args[1]
#   project_id <- args[2]
#   SAMPLE.ID <- args[3]
#   ploidy <- as.numeric(args[4])
#   cellularity <- as.numeric(args[5])
#   work_dir_cn <- args[6]
#   message("Running indels classification with the folloing settings:")
#   message("ASCAT.PATH ",ASCAT.PATH)
#   message("project_id ",project_id)
#   message("SAMPLE.ID ",SAMPLE.ID)
#   message("ploidy ",ploidy)
#   message("cellularity ",cellularity)
#   message("work_dir_cn ",work_dir_cn)
# }

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
  
  # load(file = "../lib/hrdIndex/chrominfo.RData")
  # load(file = "../lib/hrdIndex/chrominfo.snp6.hg19.RData")
  
  HRD_LOH <- calc.hrd(ascat.data2, nA=7,check.names=FALSE, return.loc=FALSE)
  # NtAI <- calc.ai(ascat.data2,chrominfo = chrominfo)[1]
  # LST <- calc.lst(ascat.data2,chrominfo = chrominfo.snp6)
  # data.frame(SampleID=SAMPLE.ID,
  #            ploidy=ploidy,
  #            cellularity=cellularity,
  #            HRD_LOH=HRD_LOH,
  #            NtAI=NtAI,
  #            LST=LST)
  # sink(file = paste0(work_dir_cn,SAMPLE.ID,".HRDscore"))
  # cat(paste0(SAMPLE.ID," ",HRD_LOH,"\n"))
  # sink()
  
  return(HRD_LOH)
}
