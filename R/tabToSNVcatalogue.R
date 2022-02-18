
#' TAB to SNV Catalogue
#' 
#' Convert a data frame containing SNV, possibly loaded from a tab seperated text file, to SNV 96 channel trinuclotide context catalogue. The data frame should containt the SNV of a single sample, and the following minimal columns: chr, position, REF, ALT.
#' 
#' @param subs data frame with subs from a single sample and the following minimal columns: chr, position, REF, ALT.
#' @param genome.v either "hg38" (will load BSgenome.Hsapiens.NCBI.GRCh38), "hg19" (will load BSgenome.Hsapiens.UCSC.hg19), mm10 (will load BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10) or canFam3 (will load BSgenome.Cfamiliaris.UCSC.canFam3::BSgenome.Cfamiliaris.UCSC.canFam3)
#' @return returns the SNV catalogue for the given sample
#' @keywords tab SNV
#' @export
#' @examples
#' subs <- read.table("subs.tab",sep="\t",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
#' res <- tabToSNVcatalogue(subs,genome.v = "hg19")
tabToSNVcatalogue <- function(subs, genome.v="hg19") {

  if(genome.v=="hg19"){
    genomeSeq <- BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
  }else if(genome.v=="hg38"){
    genomeSeq <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  }else if(genome.v=="mm10"){
    genomeSeq <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
  }else if(genome.v=="canFam3"){
    genomeSeq <- BSgenome.Cfamiliaris.UCSC.canFam3::BSgenome.Cfamiliaris.UCSC.canFam3
  }
  
  # plots mutation-context for all variants in the vcf file
  # and separately for the variants that passed
  
  mut.order <- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")
  
  #check that the required columns are present
  required_cols <- c("chr", "position", "ALT", "REF")
  if(!length(intersect(required_cols,colnames(subs)))==length(required_cols)){
    stop("[tabToSNVcatalogue error] missing columns in subs data frame, following columns required: chr, position, REF, ALT")
  }
  
  # select only SNVs
  nmutsloaded <- nrow(subs)
  selectionmuts <- nchar(as.character(subs$REF))==1 & nchar(as.character(subs$ALT))==1
  subs <- subs[selectionmuts,,drop=F]
  nmutstoanalyse <- nrow(subs)
  
  if(nmutsloaded>nmutstoanalyse){
    message("[tabToSNVcatalogue warning] mutations table contains ",nmutsloaded-nmutstoanalyse," indels, which were ignored.")
  }

  subs$wt <- subs$REF
  subs$mt <- subs$ALT
  subs$ref_base_pyrimidine_context <- subs$wt
  subs$mutant_base_pyrimidine_context <- subs$mt
  noPyrBases <- (subs$wt=='G') | (subs$wt=='A')
  subs$ref_base_pyrimidine_context[noPyrBases ] <- toPyr(subs$wt[noPyrBases])
  subs$mutant_base_pyrimidine_context[noPyrBases ]  <- toPyr(subs$mt[noPyrBases])

  # the mutations originally not in pyramidine contex
  not.pyr <- ((subs$wt=='G') | (subs$wt=='A'))
  subs$pyrwt <- as.character(subs$wt)
  subs$pyrmut <- as.character(subs$mt)
  subs$pyrwt[not.pyr] <- as.character(toPyr(subs$wt[not.pyr]))
  subs$pyrmut[not.pyr] <- as.character(toPyr(subs$mt[not.pyr]))

  # currently we support chromosomes, not contigs
  if (genome.v=="hg19"){
    expected_chroms <- paste0(c(seq(1:22),"X","Y"))
  }else if (genome.v=="hg38"){
    expected_chroms <- paste0("chr",c(seq(1:22),"X","Y"))
  }else if (genome.v=="mm10"){
    expected_chroms <- paste0("chr",c(seq(1:19),"X","Y"))
  }else if (genome.v=="canFam3"){
    expected_chroms <- paste0("chr",c(seq(1:38),"X")) 
  }
  
  # if (length(intersect(subs$chr,expected_chroms))==0) {
  #    stop("[error tabToSNVcatalogue] Input tab file does not contain seqnames ", paste(expected_chroms,collapse=" "))
  # }
  if (genome.v=="hg38" || genome.v=="mm10") {
    if(length(intersect(subs$chr,expected_chroms))==0) subs$chr <- paste0("chr",subs$chr)
  }

  muts <- data.frame(chroms=subs$chr, 
                     starts=subs$position, 
                     ends = subs$position, 
                     wt=subs$wt, 
                     mt=subs$mt, 
                     pyrwt=subs$pyrwt , 
                     pyrmut=subs$pyrmut,
                     stringsAsFactors = FALSE)


  result<- list()
  result$muts <- muts
  ######

  subs$bb <- as.character(BSgenome::getSeq(genomeSeq, as.character(subs$chr), start=subs$position-1, end=subs$position-1))
  subs$ba <- as.character(BSgenome::getSeq(genomeSeq, as.character(subs$chr), start=subs$position+1, end=subs$position+1))
  subs$wt.ref <- as.character(BSgenome::getSeq(genomeSeq, as.character(subs$chr), start=subs$position, end=subs$position))
  subs$triplets <- as.character(BSgenome::getSeq(genomeSeq, as.character(subs$chr), start=subs$position-1, end=subs$position+1))

  # table of mutations
  mut.table <- data.frame(bbef=as.character(subs$bb ), 
                          wt=as.character(subs$wt), 
                          mt=as.character(subs$mt), 
                          baft=as.character(subs$ba ),
                          stringsAsFactors = FALSE)
  mut.table$pyrwt <- mut.table$wt
  mut.table$pyrmut <- mut.table$mt
  mut.table$pyrbbef <- mut.table$bbef
  mut.table$pyrbaft <- mut.table$baft

                                  # the mutations originally not in pyramidine contex
  mut.table$pyrwt[not.pyr] <- as.character(toPyr(mut.table$wt[not.pyr]))
  mut.table$pyrmut[not.pyr] <- as.character(toPyr(mut.table$mt[not.pyr]))
  mut.table$pyrbbef[not.pyr] <- as.character(toPyr(mut.table$baft[not.pyr]))
  mut.table$pyrbaft[not.pyr] <- as.character(toPyr(mut.table$bbef[not.pyr]))

  all.hist <- generateHist(mut.table, normalise=FALSE,mut.order=mut.order)


  muts$context <- paste(mut.table$pyrbbef, '[',mut.table$pyrwt, '>',mut.table$pyrmut , ']', mut.table$pyrbaft,sep='')

  result$catalogue <- data.frame(catalogue=all.hist,row.names = names(all.hist))
  result$muts <- muts

  result

}
