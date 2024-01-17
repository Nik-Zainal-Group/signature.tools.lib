#' plot Copy Numbers from ASCAT-like file
#' 
#' Plot the copy numbers across the Chromosomes. Optionally, plot also highlight regions at the bottom. For example, could be used
#' to highlight HRD-LOH regions.
#' 
#' @param sv_df data frame with the following columns: 'seg_no', 'Chromosome', 'chromStart', 'chromEnd', 'total.copy.number.inNormal', 'minor.copy.number.inNormal', 'total.copy.number.inTumour', 'minor.copy.number.inTumour'
#' @param sample_name sample name
#' @param filename if a file name is specified, a pdf file will be generated
#' @param mar set custom margins, otherwise mar will be set to mar = c(4,4,3,2)
#' @param highlightRegions dataframe with columns 'Chromosome', 'Start', 'End', which specify the regions to highlight. This is optional.
#' @param highlightText if highlightRegions is specified, then highlightText will be plotted in the bottom right
#' @param highlightColour set the colour of the highlight regions and text. Default is brown
#' @return returns the HRD-LOH index or regions
#' @export
#' @examples
#' ascat.data <- read.table("ascat.scv",sep=",",header=TRUE)
#' HRDLOHregions <- ascatToHRDLOH(ascat.df,"test_sample",return.loc=TRUE)
#' HRDLOHscore <- nrow(HRDLOHregions)
#' plotCopyNumbers(sv_df = ascat.data,sample_name = "",highlightRegions = HRDLOHregions,
#'                 highlightText = paste0("HRD-LOH\\n",HRDLOHscore))
plotCopyNumbers <- function(sv_df,
                            sample_name,
                            filename=NULL,
                            plottitle=NULL,
                            mar=NULL,
                            minorCNcolour="green",
                            totalCNcolour="red",
                            outofrangeCNcolour="purple",
                            highlightRegions=NULL,
                            highlightText=NULL,
                            highlightColour = "brown",
                            textscaling = 1,
                            genome.v="hg19"){

  if(!is.null(filename)) cairo_pdf(filename = filename,width = 7,height = 3)
  if(!is.null(mar)) {
    par(mar=mar)
  } else {
    mar = c(4,4,3,2)
  }
  if(is.null(plottitle)) plottitle <- sample_name
  if(genome.v=="hg38"){
    chromosome_lengths <- GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[1:24]
    chrom_names <- substr(names(chromosome_lengths),4,5)
  }else if(genome.v=="hg19"){
    chromosome_lengths <- GenomeInfoDb::seqlengths(BSgenome.Hsapiens.1000genomes.hs37d5::hs37d5)[1:24]
    chrom_names <- names(chromosome_lengths)
  }else if(genome.v=="mm10"){
    chromosome_lengths <- GenomeInfoDb::seqlengths(BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)[1:21]
    chrom_names <- substr(names(chromosome_lengths),4,5)
  }else if(genome.v=="canFam3"){
    chromosome_lengths <- GenomeInfoDb::seqlengths(BSgenome.Cfamiliaris.UCSC.canFam3::BSgenome.Cfamiliaris.UCSC.canFam3)[1:39]
    chrom_names <- substr(names(chromosome_lengths),4,5)
  }else{
    message("[error plotCopyNumbers] genome version ",genome.v," not supported.")
    return(NULL)
  }

  displace <- 0.05
  barheight <- 0.1
  highestCN <- min(max(sv_df$total.copy.number.inTumour)+1,6) 
  maxCoord <- sum(as.numeric(chromosome_lengths))
  plot(NA,xlim=c(0,maxCoord),ylim=c(0-0.3,highestCN+0.6),main=plottitle,cex.main=0.8,
       ylab = "Copy Numbers",xlab = "",xaxt="n",xaxs="i",yaxs="i",mgp=c(1.5,0.6,0),mar=mar,
       cex = textscaling,cex.lab = textscaling, cex.main = textscaling,cex.axis = textscaling,las = 2)
  # abline(v=maxCoord)
  for(i in 0:floor(highestCN)) lines(x=c(0,maxCoord),y=c(i,i),col="lightgrey",lty=3,lwd=0.5)
  for (i in 1:length(chromosome_lengths)){
    chrom_name <- chrom_names[i]
    tmpChrom <- sv_df[sv_df$Chromosome==chrom_name,]
    startCoord <- ifelse(i==1,0,sum(as.numeric(chromosome_lengths[1:(i-1)])))
    if(startCoord>0) abline(v=startCoord)
    par(xpd=TRUE)
    text(x=startCoord+(chromosome_lengths[i])/2,y=-1.2*textscaling,labels = chrom_name,cex = 0.6*textscaling)
    par(xpd=FALSE)
    if(nrow(tmpChrom)>0){
      for (j in 1:nrow(tmpChrom)){
        if(tmpChrom$total.copy.number.inTumour[j]<7){
          rect(xleft = startCoord+tmpChrom$chromStart[j],xright = startCoord+tmpChrom$chromEnd[j],
               ybottom = tmpChrom$total.copy.number.inTumour[j]+displace,
               ytop = tmpChrom$total.copy.number.inTumour[j]+displace+barheight,col = totalCNcolour,border = totalCNcolour)
        }else{
          rect(xleft = startCoord+tmpChrom$chromStart[j],xright = startCoord+tmpChrom$chromEnd[j],
               ybottom = 6 + 4*displace+barheight,
               ytop = 6 + 4*displace+2*barheight,col = outofrangeCNcolour,border = outofrangeCNcolour)
        }     
        rect(xleft = startCoord+tmpChrom$chromStart[j],xright = startCoord+tmpChrom$chromEnd[j],
             ybottom = tmpChrom$minor.copy.number.inTumour[j]-displace-barheight,
             ytop = tmpChrom$minor.copy.number.inTumour[j]-displace,col = minorCNcolour,border = minorCNcolour)
      }
    }
    # add highlight regions
    if(!is.null(highlightRegions)){
      tmpChromH <- highlightRegions[highlightRegions$Chromosome==chrom_name,]
      if(nrow(tmpChromH)>0){
        for (j in 1:nrow(tmpChromH)){
          par(xpd=TRUE)
          rect(xleft = startCoord + tmpChromH$Start[j],xright = startCoord + tmpChromH$End[j],
               ytop = -0.15,ybottom = -0.75,col = highlightColour,border = highlightColour)
          par(xpd=FALSE)
        }
      }
    }
  }
  par(xpd=TRUE)
  if(!is.null(highlightRegions)){
    rect(xleft = 0,xright = maxCoord,ytop = -0.15,ybottom = -0.75)
    if(!is.null(highlightText)) text(x=2*-10^8,y=-0.05,labels = highlightText,cex = 0.8*textscaling,col = highlightColour,pos = 1)
  }
  text(x=maxCoord/2,y=-2*textscaling,labels = "chromosomes",cex = textscaling)
  par(xpd=FALSE)
  legend("topleft",legend = c("Minor","Total","Total>6"),fill = c(minorCNcolour,totalCNcolour,outofrangeCNcolour),border = rep("white",3),bty = "n",horiz = TRUE,xpd = TRUE,inset = c(0,-0.3),cex=textscaling*0.8)
  if(!is.null(filename)) dev.off()

}                                                                                                           
                                                                                                             