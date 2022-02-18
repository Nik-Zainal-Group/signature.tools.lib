
#' VCF to SNV catalogue
#' 
#' Convert a vcf file containing SNV to SNV 96 channel trinuclotide context catalogue. The VCF file should containt the SNV of a single sample.
#' 
#' @param vcfFilename path to input VCF (file must be tabix indexed)
#' @param genome.v either "hg38" (will load BSgenome.Hsapiens.UCSC.hg38), "hg19" (will load BSgenome.Hsapiens.1000genomes.hs37d5), mm10 (will load BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10) or canFam3 (will load BSgenome.Cfamiliaris.UCSC.canFam3::BSgenome.Cfamiliaris.UCSC.canFam3)
#' @return returns the SNV catalogue for the given sample
#' @keywords vcf SNV
#' @export
#' @examples
#' file_subs <- "subs.vcf"
#' res <- vcfToSNVcatalogue(file_subs,genome.v = "hg38")
vcfToSNVcatalogue <- function(vcfFilename, genome.v="hg19") {
  
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
 
  mut.order <- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")
  
  # read only chr seqnames from VCF, not contigs
  #gr <- GenomicRanges::GRanges(GenomeInfoDb::Seqinfo(genome=genome.v))
  gr <- GenomicRanges::GRanges(GenomeInfoDb::seqinfo(genomeSeq))
  # if (genome.v=="hg38" || genome.v=="mm10") {
  #   GenomeInfoDb::seqlevels(gr) <- sub("chr", "", GenomeInfoDb::seqlevels(gr))
  # }
  vcf_seqnames <- Rsamtools::headerTabix(vcfFilename)$seqnames 
    if (genome.v=="hg38" || genome.v=="mm10") {
    if(length(intersect(vcf_seqnames,expected_chroms))==0) vcf_seqnames <- paste0("chr",vcf_seqnames)
  }
  
  if(tools:::.BioC_version_associated_with_R_version()<3.5){
    gr <- GenomeInfoDb::keepSeqlevels(gr,intersect(vcf_seqnames,expected_chroms))
  }else{
    gr <- GenomeInfoDb::keepSeqlevels(gr,intersect(vcf_seqnames,expected_chroms),pruning.mode = "coarse")
  }
  
  vcf_seqnames <- Rsamtools::headerTabix(vcfFilename)$seqnames
  if (genome.v=="hg38" || genome.v=="mm10") {
    if(length(intersect(vcf_seqnames,expected_chroms))==0) GenomeInfoDb::seqlevels(gr) <- sub("chr", "", GenomeInfoDb::seqlevels(gr))
  }
  
  # load the vcf file
  vcf_data <- VariantAnnotation::readVcf(vcfFilename, genome.v, gr)
  vcf_data <- VariantAnnotation::expand(vcf_data)
  nmutsloaded <- nrow(vcf_data)
  
  #filters failed for each variant
  rd <- SummarizedExperiment::rowRanges(vcf_data)
  selectionvcf <- nchar(as.character(rd$REF))==1 & nchar(as.character(rd$ALT))==1
  vcf_data <- vcf_data[selectionvcf,,drop=F]
  nmutssnv <- nrow(vcf_data)
  
  if(nmutsloaded>nmutssnv){
    message("[vcfToSNVcatalogue warning] the vcf file ",vcfFilename," contains ",nmutsloaded-nmutssnv," indels, which were ignored.")
  }
  
  rd <- SummarizedExperiment::rowRanges(vcf_data)
  info.data <- VariantAnnotation::info(vcf_data)
  
  rgs <- IRanges::ranges(vcf_data)
  starts <- BiocGenerics::start(rgs)
  ends <-  BiocGenerics::end(rgs)
  
  #Check chromosomes exist
  chroms <- GenomeInfoDb::seqnames(vcf_data)

  if (length(chroms)==0){ 
     stop("[vcfToSNVcatalogue error] Input vcf does not contain variants ", vcfFilename)
  }
  
  if (genome.v=="hg38" || genome.v=="mm10") {
    if(length(intersect(chroms,expected_chroms))==0) chroms <- paste0("chr",chroms)
  }
  
  fxd <- (VariantAnnotation::fixed(vcf_data))
  wt <- as.character(rd$REF)
  mt <- as.character(rd$ALT)
  
  
  barcode <- paste(chroms, '-',starts,'-', mt, sep='')
  
  #get context
  bb <- as.character(BSgenome::getSeq(genomeSeq, chroms, start=starts-1, end=ends-1))
  ba <- as.character(BSgenome::getSeq(genomeSeq, chroms, start=starts+1, end=ends+1))
  wt.ref <- as.character(BSgenome::getSeq(genomeSeq, chroms, start=starts, end=ends))
  triplets <- as.character(BSgenome::getSeq(genomeSeq, chroms, start=starts-1, end=ends+1))

  mut.table <- data.frame(bbef=as.character(bb), wt=as.character(wt), mt=as.character(mt), baft=as.character(ba), stringsAsFactors=FALSE)


  mut.table$pyrwt <- as.character(mut.table$wt)
  mut.table$pyrmut <- as.character(mut.table$mt)
  mut.table$pyrbbef <- as.character(mut.table$bbef)
  mut.table$pyrbaft <- as.character(mut.table$baft)


  # the mutations originally not in pyramidine contex
  not.pyr <- ((wt=='G') | (wt=='A'))
  mut.table$pyrwt[not.pyr] <- as.character(toPyr(mut.table$wt[not.pyr]))
  mut.table$pyrmut[not.pyr] <- as.character(toPyr(mut.table$mt[not.pyr]))
  mut.table$pyrbbef[not.pyr] <- as.character(toPyr(mut.table$baft[not.pyr]))
  mut.table$pyrbaft[not.pyr] <- as.character(toPyr(mut.table$bbef[not.pyr]))

  all.hist <- generateHist(mut.table, normalise=FALSE,mut.order=mut.order)

  # names(all.hist) <- mut.order

  muts <- data.frame(chroms=chroms, starts=starts, ends = ends, wt=wt, mt=mt, pyrwt=mut.table$pyrwt , pyrmut=mut.table$pyrmut, pass=TRUE, barcode=barcode,
                 context=paste(mut.table$pyrbbef, '[',mut.table$pyrwt, '>',mut.table$pyrmut , ']', mut.table$pyrbaft,sep=''),
                 tumor.freq=rep(NA,length(chroms)), normal.freq=rep(NA,length(chroms)),
                 tumor.reads=rep(NA,length(chroms)), normal.reads=rep(NA,length(chroms)),
                 tumor.depth=rep(NA,length(chroms)), normal.depth=rep(NA,length(chroms)),
                                                         isSnp =rep(NA,length(chroms)),
                 stringsAsFactors = F)
  
  result<- list()
  result$catalogue <- data.frame(catalogue=all.hist,row.names = names(all.hist))
  result$muts <- muts
  result$mut.table <- mut.table
  result

}
