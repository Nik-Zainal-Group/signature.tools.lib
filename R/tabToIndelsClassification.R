
#' tab to Indels Classification
#' 
#' Convert a dataframe containing Indels into a data frame where each indel is classified as repet-mediated, Microhomology-mediated or other. A summary of the count of indels (deletions and insertions) and their proportion with respect to the total is also provided.
#' 
#' @param indel.data dataframe with indels from a single sample and the following minimal columns: chr, position, REF, ALT.  
#' @param sampleID name of the sample
#' @param genome.v version of the genome to be used to look up the context of the indel, either "hg19", "hg38", "mm10" or "canFam3"
#' @return the function returns a list with elements "indels_classified", which is a table with the indels and their classification, and "count_proportion", which is a summary of the count of indels and their proportion
#' @export
#' @examples 
#' res <- tabToIndelsClassification(indel.data,"testSample","hg19")
tabToIndelsClassification <- function(indel.data,
                                      sampleID,
                                      genome.v="hg19"){
  
  if(genome.v=="hg19"){
    expected_chroms <- paste0(c(seq(1:22),"X","Y"))
    genomeSeq <- BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
  }else if(genome.v=="hg38"){
    expected_chroms <- paste0("chr",c(seq(1:22),"X","Y"))
    genomeSeq <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  }else if(genome.v=="mm10"){
    expected_chroms <- paste0("chr",c(seq(1:19),"X","Y"))
    genomeSeq <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
  }else if(genome.v=="canFam3"){
    expected_chroms <- paste0("chr",c(seq(1:38),"X")) 
    genomeSeq <- BSgenome.Cfamiliaris.UCSC.canFam3::BSgenome.Cfamiliaris.UCSC.canFam3
  }
  
  # read only chr seqnames from VCF, not contigs
  vcf_seqnames <- unique(indel.data$chr) 
  if (genome.v=="hg38" | genome.v=="mm10" | genome.v=="canFam3") {
    if(length(intersect(vcf_seqnames,expected_chroms))==0) indel.data$chr <- paste0("chr",indel.data$chr)
  }
  indel.data <- indel.data[indel.data$chr %in% expected_chroms,]
  
  if(nrow(indel.data)==0){
    message("[warning tabToIndelsClassification] no indels founds, nothing to process.")
    return(NULL)
  }
  
  # convert formats, and find context of the indels
  indel.df <- prepare.indel.df_tabversion(indel.data,genomeSeq,genome.v)
  
  indels_classified <- mh(indel.df)
  
  # let's add another column to clarify the class of each indel
  indels_classified$indel.class <- "-"
  indels_classified$indel.class[indels_classified$classification=="Microhomology-mediated"] <- "del.mhomology"
  indels_classified$indel.class[indels_classified$classification=="Repeat-mediated"] <- "del.repeatmediated"
  indels_classified$indel.class[indels_classified$classification=="None"] <- "del.other"
  indels_classified$indel.class[indels_classified$indel.type=="I"] <- "insertion"
  indels_classified$indel.class[indels_classified$indel.type=="DI"] <- "indel.complex"
  
  # save and return
  res <- list()
  res$indels_classified <- indels_classified
  res$count_proportion <- indelsToCountAndProportion(res$indels_classified,sampleID)
  
  return(res)
  
}

###########################################################

prepare.indel.df_tabversion <- function(indel.data,genomeSeq,genome.v) {
  
  if (nrow(indel.data)>0) {
    
    ref.length <- nchar(indel.data$REF)
    alt.length <- nchar(indel.data$ALT)
    indel.length <- abs(ref.length - alt.length)
    
    indel.type <- rep(NA, nrow(indel.data))
    indel.type[ref.length==1 & alt.length>1] <- 'I'
    indel.type[ref.length>1 & alt.length==1] <- 'D'
    indel.type[ref.length>1 & alt.length>1] <- 'DI'
    indel.type[ref.length==1 & alt.length==1] <- 'DI'
    
    # sequence of change
    change <- vector()
    change[indel.type=='DI'] <-  substr( as.character(indel.data$REF)[indel.type=='DI'],2,1e5)
    change[indel.type=='I'] <- substr( as.character(indel.data$ALT)[indel.type=='I'], 2, 1e5)
    change[indel.type=='D'] <- substr( as.character(indel.data$REF), 2, 1e5)[indel.type=='D']
    
    min.position <- indel.data$position
    max.position <- indel.data$position + indel.length 
    indel.chr <- as.character(indel.data$chr)

    extend5 = min.position-indel.length-25;
    extend3 = max.position + indel.length+25;
    
    
    slice5 <- as.character(BSgenome::getSeq(genomeSeq, indel.chr, extend5, min.position))
    # 
    slice3 <- as.character(BSgenome::getSeq(genomeSeq, indel.chr, max.position+1, extend3))
    
    indel.df <- data.frame(
      chr=as.character(indel.data$chr),
      pos=indel.data$position,
      ref=as.character(indel.data$REF),
      alt=as.character(indel.data$ALT),
      indel.type=indel.type,
      change=change,
      slice3=slice3,
      slice5=slice5,
      indel.length=indel.length
    ) } else {
      
      indel.df <- data.frame()
      
      
    }
  
  indel.df
}
