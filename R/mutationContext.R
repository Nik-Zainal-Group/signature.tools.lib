
#' Find SNV context
#'
#' Given a table of SNVs, find the reference genome nucleotide bases that are at
#' the 5' and 3' of each SNV location, according to the pyrimidine reference. For
#' example, if the mutation is a TT[C>T]AG, then the 5' will correspond to lower genomic
#' coordinates, with 5' context TT and the 3' will correspond to the higher genomic coordinates,
#' with 3' context AG. Alternatively, if the mutation is a TT[A>T]AG, then the reference
#' considered will be on the opposite strand, and the mutation and contexts will be
#' reversed as CT[T>A]AA.
#' It is also possible to filter the mutations in the given SNV table according to
#' the specific pyrimidine mutation and, optionally, the first base on the 5' or 3'.
#' 
#' 
#' @param snv_table data frame containing SNVs, with required columns chr, position, REF, ALT
#' @param mtype if specified, filter mutations according to specific base changes: C>A, C>T, C>G, T>A, T>C, T>G 
#' @param fiveprime if specified, filter mutations according to specific 5': C, A, T, G. Used only if mtype is specified
#' @param fiveprime if specified, filter mutations according to specific 3': C, A, T, G. Used only if mtype is specified
#' @param context_length number of nucleotides on the 5' and 3' to consider
#' @param genomev genome version, hg19 or hg38
#' @return table of counts of base contexts
#' @export
findContextSNV <- function(snv_table,
                           mtype = NULL,
                           fiveprime = NULL,
                           threeprime = NULL,
                           genomev = "hg19",
                           context_length = 5){
  #now use the catalogue function to get the trinucletide context
  cat_res <- tabToSNVcatalogue(subs = snv_table,genome.v = genomev)
  annotated_SNV <- cat_res$muts
  #not get me only the N[T>G]G mutations
  if(!is.null(mtype)){
    mtype <- paste0("[",mtype,"]")
    startp <- 2
    stopp <- 6
    if(!is.null(fiveprime)){
      mtype <- paste0(fiveprime,mtype)
      startp <- startp - 1
    }
    if(!is.null(threeprime)){
      mtype <- paste0(mtype,threeprime)
      stopp <- stopp + 1
    }
    select_muts <-  substr(annotated_SNV$context,start = startp,stop = stopp) == mtype
    annotated_SNV <- annotated_SNV[select_muts,]
  }

  #load the genome
  if (genomev == "hg19") {
    genomeSeq <- BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
  }else if (genomev == "hg38") {
    genomeSeq <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  }
  #tell me the context
  howmany <- context_length
  five_prime <- as.character(BSgenome::getSeq(genomeSeq, as.character(annotated_SNV$chroms), start=annotated_SNV$starts-howmany, end=annotated_SNV$ends-1))
  three_prime <- as.character(BSgenome::getSeq(genomeSeq, as.character(annotated_SNV$chroms), start=annotated_SNV$starts+1, end=annotated_SNV$ends+howmany))
  
  #now let's count and plot nicely
  five_prime_counts <- matrix(0,nrow = howmany,ncol = 4,
                              dimnames = list(1:howmany,c("C","T","A","G")))
  three_prime_counts <- matrix(0,nrow = howmany,ncol = 4,
                               dimnames = list(1:howmany,c("C","T","A","G")))
  
  flip <- function(x){
    char_list <- strsplit(x,split="")[[1]]
    return_str <- c()
    for (i in 1:length(char_list)) {
      if(char_list[i]=="A"){
        return_str <- c("T",return_str)
      }else if(char_list[i]=="C"){
        return_str <- c("G",return_str)
      }else if(char_list[i]=="T"){
        return_str <- c("A",return_str)
      }else if(char_list[i]=="G"){
        return_str <- c("C",return_str)
      } 
    }
    paste0(return_str,collapse = "")
  }
  
  #which need to be flipped?
  to_be_flipped <- annotated_SNV$wt != annotated_SNV$pyrwt
  for (i in 1:nrow(annotated_SNV)) {
    if (to_be_flipped[i]){
      tmp <- flip(five_prime[i])
      five_prime[i] <- flip(three_prime[i])
      three_prime[i] <- tmp
    }
  }
  
  for (i in 1:howmany){
    for (r in 1:nrow(annotated_SNV)){
      b <- substr(five_prime[r],i,i)
      five_prime_counts[i,b] <- five_prime_counts[i,b] + 1
      b <- substr(three_prime[r],i,i)
      three_prime_counts[i,b] <- three_prime_counts[i,b] + 1
    }
  }
  
  #combine
  pyrRefcount <- table(annotated_SNV$pyrwt)
  t_row <- c(0,0,0,0)
  names(t_row) <- c("C","T","A","G")
  t_row[names(pyrRefcount)] <- as.vector(pyrRefcount)
  data_table <- rbind(as.data.frame(five_prime_counts),t_row,as.data.frame(three_prime_counts))
  return(data_table)
}

#' plot SNV context
#'
#' Plot the result table obtained from calling the findContextSNV function
#' 
#' 
#' @param data_table table with counts of nucleotides at the 5' and 3' of a mutation, obtained from calling the findContextSNV function
#' @param outfile if specified the plot will be save to a pdf file, use extension .pdf
#' @param main optional plot title
#' @param xlab optional x label
#' @export
plotContextSNV <- function(data_table,
                           outfile = NULL,
                           main="",
                           xlab=""){
  howmany <- floor(nrow(data_table)/2)
  # plot
  kelly_colors <- c('F3C300', '875692', 'F38400', 'A1CAF1', 'BE0032', 
                    'C2B280', '848482', '008856', 'E68FAC', '0067A5', 'F99379', '604E97', 
                    'F6A600', 'B3446C', 'DCD300', '882D17', '8DB600', '654522', 'E25822', '2B3D26','222222','F2F3F4', 'CCCCCC','CCCCCC','CCCCCC')
  kelly_colors <- paste0("#",kelly_colors)
  if(!is.null(outfile)) cairo_pdf(outfile,width = 10,height = 6)
  par(mar=c(5,5,3,4))
  barplot(t(as.matrix(data_table)),names.arg = -howmany:howmany,col = kelly_colors[c(3,4,1,2)],
          main =  main,
          ylab = "count",
          xlab = xlab)
  par(xpd=TRUE)
  legend("topright",inset=c(-0.05,0),legend = c("C","T","A","G"),fill = kelly_colors[c(3,4,1,2)],border = NA,bty = "n")
  par(xpd=FALSE)
  abline(h = sum(data_table[1,])/2,lty = 2)
  if(!is.null(outfile)) dev.off()
}

