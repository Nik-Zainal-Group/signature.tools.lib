
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
tabToIndelsClassification <- function(indel.data,sampleID, genome.v="hg19"){
  
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
  #gr <- GenomicRanges::GRanges(GenomeInfoDb::Seqinfo(genome=genome.v))
  gr <- GenomicRanges::GRanges(GenomeInfoDb::seqinfo(genomeSeq))
  if (genome.v=="hg38" || genome.v=="mm10") {
    GenomeInfoDb::seqlevels(gr) <- sub("chr", "", GenomeInfoDb::seqlevels(gr))
  }
  vcf_seqnames <- unique(indel.data$chr) 
  
  if(tools:::.BioC_version_associated_with_R_version()<3.5){
    gr <- GenomeInfoDb::keepSeqlevels(gr,intersect(vcf_seqnames,expected_chroms))
  }else{
    gr <- GenomeInfoDb::keepSeqlevels(gr,intersect(vcf_seqnames,expected_chroms),pruning.mode = "coarse")
  }
  
  # convert formats, and find context of the indels
  indel.df <- prepare.indel.df_tabversion(indel.data,genomeSeq,genome.v,expected_chrom)
  
  res <- list()
  res$indels_classified <- mh(indel.df)
  res$count_proportion <- indelsToCountAndProportion(res$indels_classified,sampleID)
  
  return(res)
  
}

###########################################################

prepare.indel.df_tabversion <- function(indel.data,genomeSeq,genome.v,expected_chrom) {
  
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
    # if (genomeSeq@provider_version=="mm10"){
    #   indel.chr <- paste('chr',indel.chr,sep='')
    # }
    if (genome.v=="hg38" || genome.v=="mm10") {
      if(length(intersect(indel.chr,expected_chroms))==0) indel.chr <- paste0("chr",indel.chr)
    }
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
