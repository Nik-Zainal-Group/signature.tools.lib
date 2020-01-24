
#' VCF to Indels Classification
#' 
#' Convert a VCF file containing Indels into a data frame where each indel is classified as repet-mediated, Microhomology-mediated or other. A summary of the count of indels (deletions and insertions) and their proportion with respect to the total is also provided.
#' 
#' @param indelsVCF.file path to input VCF (file must be tabix indexed). This file should have been already filtered for the final indels sets to be used in the analysis. 
#' @param sampleID name of the sample
#' @param genome.v version of the genome to be used to look up the context of the indel, either "hg19", "hg38", "mm10" or "canFam3"
#' @return the function returns a list with elements "indels_classified", which is a table with the indels and their classification, and "count_proportion", which is a summary of the count of indels and their proportion
#' @export
#' @examples 
#' res <- vcfToIndelsClassification("test.indel.vcf.gz","testSample","hg19")
vcfToIndelsClassification <- function(indelsVCF.file,sampleID, genome.v="hg19"){

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
  # if (genome.v=="hg38" || genome.v=="mm10") {
  #   GenomeInfoDb::seqlevels(gr) <- sub("chr", "", GenomeInfoDb::seqlevels(gr))
  # }
  vcf_seqnames <- Rsamtools::headerTabix(indelsVCF.file)$seqnames 
  if (genome.v=="hg38" || genome.v=="mm10") {
    if(length(intersect(vcf_seqnames,expected_chroms))==0) vcf_seqnames <- paste0("chr",vcf_seqnames)
  }
  
  if(tools:::.BioC_version_associated_with_R_version()<3.5){
    gr <- GenomeInfoDb::keepSeqlevels(gr,intersect(vcf_seqnames,expected_chroms))
  }else{
    gr <- GenomeInfoDb::keepSeqlevels(gr,intersect(vcf_seqnames,expected_chroms),pruning.mode = "coarse")
  }
    
  vcf_seqnames <- Rsamtools::headerTabix(indelsVCF.file)$seqnames
  if (genome.v=="hg38" || genome.v=="mm10") {
    if(length(intersect(vcf_seqnames,expected_chroms))==0) GenomeInfoDb::seqlevels(gr) <- sub("chr", "", GenomeInfoDb::seqlevels(gr))
  }
  
  # load the indel VCF file
  indel.data <- VariantAnnotation::readVcf(indelsVCF.file, genome.v, gr)
  indel.data <- VariantAnnotation::expand(indel.data)
  
  # convert formats, and find context of the indels
  indel.df <- prepare.indel.df(indel.data,genomeSeq,genome.v,expected_chroms)
  
  res <- list()
  res$indels_classified <- mh(indel.df)
  res$count_proportion <- indelsToCountAndProportion(res$indels_classified,sampleID)
  
  return(res)

}

###########################################################

prepare.indel.df <- function(indel.data,genomeSeq,genome.v,expected_chroms) {
  
  if (nrow(indel.data)>0) {
    
    
    ref.length <- Biostrings::width(SummarizedExperiment::rowRanges(indel.data)$REF)
    alt.length <- Biostrings::width(SummarizedExperiment::rowRanges(indel.data)$ALT)
    indel.length <- abs(ref.length - alt.length)
    
    indel.type <- rep(NA, nrow(indel.data))
    indel.type[ref.length==1 & alt.length>1] <- 'I'
    indel.type[ref.length>1 & alt.length==1] <- 'D'
    indel.type[ref.length>1 & alt.length>1] <- 'DI'
    indel.type[ref.length==1 & alt.length==1] <- 'DI'
    
    # sequence of change
    change <- vector()
    change[indel.type=='DI'] <-  substr( as.character(SummarizedExperiment::rowRanges(indel.data)$REF)[indel.type=='DI'],2,1e5)
    change[indel.type=='I'] <- substr( as.character(SummarizedExperiment::rowRanges(indel.data)$ALT)[indel.type=='I'], 2, 1e5)
    change[indel.type=='D'] <- substr( as.character(SummarizedExperiment::rowRanges(indel.data)$REF), 2, 1e5)[indel.type=='D']
    
    min.position <- BiocGenerics::start(indel.data)
    max.position <- BiocGenerics::start(indel.data) + indel.length 
    indel.chr <- as.character(GenomeInfoDb::seqnames(indel.data))
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
      chr=as.character(GenomeInfoDb::seqnames(indel.data)),
      pos=BiocGenerics::start(IRanges::ranges(indel.data)),
      ref=as.character(SummarizedExperiment::rowRanges(indel.data)$REF),
      alt=as.character(SummarizedExperiment::rowRanges(indel.data)$ALT),
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


##############################################

indelsToCountAndProportion <- function(all.indels.table, sampleIDs) {
  
  if(nrow(all.indels.table)>0){
    if (length(sampleIDs)==1) {
      all.indels.table$sample <- sampleIDs
    }
    
    
    all.deletions.table.mh <- subset(all.indels.table, indel.type=='D' & classification=='Microhomology-mediated')
    if (nrow(all.deletions.table.mh)>0) {
      deletions.mh.samples <- as.data.frame(table(all.deletions.table.mh$sample))
      colnames(deletions.mh.samples) <- c('sample', 'del.mh')
      rownames(deletions.mh.samples) <- deletions.mh.samples$sample
    } else {
      deletions.mh.samples <- data.frame(sample=sampleIDs, del.mh=0)
      rownames(deletions.mh.samples) <- sampleIDs
    }
    
    deletions.table.repeat <- subset(all.indels.table, indel.type=='D' & classification=='Repeat-mediated')
    if (nrow(deletions.table.repeat)>0) {
      deletions.repeat.samples <- as.data.frame(table(deletions.table.repeat$sample))
      colnames(deletions.repeat.samples) <- c('sample', 'del.rep')
      rownames(deletions.repeat.samples) <- deletions.repeat.samples$sample
    } else {
      deletions.repeat.samples <- data.frame(sample=sampleIDs, del.rep=0)
      rownames(deletions.repeat.samples) <- sampleIDs
    }
    
    all.deletions.table.other <- subset(all.indels.table, indel.type=='D' & classification=='None')
    if (nrow(all.deletions.table.other)) {
      deletions.other.samples <- as.data.frame(table(all.deletions.table.other$sample))
      colnames(deletions.other.samples ) <- c('sample', 'del.other')
      rownames(deletions.other.samples) <- deletions.other.samples$sample
    } else {
      deletions.other.samples <- data.frame(sample=sampleIDs, del.other=0)
      rownames(deletions.other.samples) <- sampleIDs
    }
    
    all.insertions.table <- subset(all.indels.table, indel.type=='I')
    if (nrow(all.insertions.table) >0) {
      insertions.samples <- as.data.frame(table(all.insertions.table $sample))
      colnames(insertions.samples ) <- c('sample', 'ins')
      rownames(insertions.samples) <- insertions.samples$sample
    } else {
      insertions.samples <- data.frame(sample=sampleIDs, ins=0)
      rownames(insertions.samples) <- sampleIDs
    }
    
    indel.table <- data.frame(sample=sampleIDs,
                              del.mh=deletions.mh.samples[as.character(sampleIDs), 'del.mh'],
                              del.rep=deletions.repeat.samples[as.character(sampleIDs), 'del.rep'],
                              del.none=deletions.other.samples[as.character(sampleIDs), 'del.other'],
                              ins=insertions.samples[as.character(sampleIDs), 'ins']
    )
    
    indel.table$all.indels <- indel.table$del.mh+indel.table$del.rep+indel.table$del.none
    
    
    # del.mh.prop del.rep.prop del.none.prop
    indel.table$del.mh.prop <- indel.table$del.mh/indel.table$all.indels
    indel.table$del.rep.prop <- indel.table$del.rep/indel.table$all.indels
    indel.table$del.none.prop <- indel.table$del.none/indel.table$all.indels
    
    indel.table$del.mh.count <- indel.table$del.mh
    indel.table$del.rep.count <- indel.table$del.rep
    indel.table$del.none.count <- indel.table$del.none
  }else{
    indel.table <- data.frame(sample=sampleIDs,
                              del.mh=rep(0,length(sampleIDs)),
                              del.rep=rep(0,length(sampleIDs)),
                              del.none=rep(0,length(sampleIDs)),
                              ins=rep(0,length(sampleIDs)))
    
    indel.table$all.indels <- indel.table$del.mh+indel.table$del.rep+indel.table$del.none
    
    
    # del.mh.prop del.rep.prop del.none.prop
    indel.table$del.mh.prop <- indel.table$del.mh
    indel.table$del.rep.prop <- indel.table$del.rep
    indel.table$del.none.prop <- indel.table$del.none
    
    indel.table$del.mh.count <- indel.table$del.mh
    indel.table$del.rep.count <- indel.table$del.rep
    indel.table$del.none.count <- indel.table$del.none
  }
  
  indel.table
  
}

