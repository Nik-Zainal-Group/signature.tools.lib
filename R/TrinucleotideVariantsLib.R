
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

#' Build a Triinucleotide Variants Catalogue from DNVs
#' 
#' This function takes as input a list of double nucleotide variants (DNVs),
#' and computes a list of trinucleotide variants (TNVs) finding which DNVs
#' are next to each other. It then returns the annotated TNVs and the TNV catalogues.
#' The trinucleotide mutations considered are only those where all three bases change,
#' i.e. TTT>AAA but not TTT>ATA.
#' 
#' @param dnvtab requires columns Sample, Chrom, Pos, Ref, Alt, with Ref and Alt of length 2
#' @return list of TNVs and TNV catalogue
#' @export
dnvTabToTNVcatalogue <- function(dnvtab){
  
  # channels
  bases <- c("A","C","G","T")
  combos <- c()
  for (b1 in bases){
    for (b2 in bases){
      for (b3 in bases){
        combos <- c(combos,paste0(b1,b2,b3))
      }
    }
  }
  combos_rc <- sapply(combos,flip,USE.NAMES = F)
  # choose complement by alphabetical order
  chosen_trinucl <- unique(sapply(1:length(combos),function(i)ifelse(combos[i]>combos_rc[i],combos[i],combos_rc[i]),USE.NAMES = F))
  # now to get all mutation types
  mutationTypes <- c()
  for (trinucl in chosen_trinucl) {
    t1 <- substr(trinucl,1,1)
    t2 <- substr(trinucl,2,2)
    t3 <- substr(trinucl,3,3)
    for (possibleMut in combos) {
      # check if any base is the same
      m1 <- substr(possibleMut,1,1)
      m2 <- substr(possibleMut,2,2)
      m3 <- substr(possibleMut,3,3)
      if(t1 != m1 & t2 != m2 & t3 != m3) mutationTypes <- c(mutationTypes,paste(trinucl,possibleMut,sep = ">"))
    }
  }
  
  sample_list <- unique(dnvtab$Sample)
  TNV_catalogue <- data.frame(row.names = mutationTypes,stringsAsFactors = F)
  TNV_table <- NULL
  for (s in sample_list){
    # for debug: s <- sample_list[1]
    sample_dnvs <- dnvtab[dnvtab$Sample==s,,drop=F]
    allchroms <- unique(sample_dnvs$Chrom)
    sample_TNVs <- NULL
    for (chrom in allchroms){
      # for debug: chrom <- allchroms[1]
      # search for DNV in each chromosome
      chromsubs <- sample_dnvs[sample_dnvs$Chrom==chrom,,drop=F]
      # order by position
      chromsubs <- chromsubs[order(chromsubs$Pos),,drop=F]
      # now annotate
      chromsubs$Pos_nextDNV <- c(chromsubs$Pos[-1],0)
      chromsubs$dist_nextDNV <- chromsubs$Pos-chromsubs$Pos_nextDNV
      chromsubs$Ref_nextDNV <- c(chromsubs$Ref[-1],"NN")
      chromsubs$Alt_nextDNV <- c(chromsubs$Alt[-1],"NN")
      # now find the TNVs
      trinuc_index <- which(chromsubs$dist_nextDNV==-1)
      chromsubs <- chromsubs[trinuc_index,,drop=F]
      if (nrow(chromsubs)>0){
        # if some TNVs are found, annotate
        chromsubs$trinuc_Ref <- sapply(1:nrow(chromsubs),function(i) paste0(chromsubs$Ref[i],substr(chromsubs$Ref_nextDNV[i],2,2)))
        chromsubs$trinuc_Alt <- sapply(1:nrow(chromsubs),function(i) paste0(chromsubs$Alt[i],substr(chromsubs$Alt_nextDNV[i],2,2)))
        chromsubs$trinuc_mutation <- paste0(chromsubs$trinuc_Ref,">",chromsubs$trinuc_Alt)
        # get reverse complement
        # chromsubs$trinuc_Ref_rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(chromsubs$dinuc_Ref)))
        # chromsubs$trinuc_Alt_rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(chromsubs$dinuc_Alt)))
        chromsubs$trinuc_Ref_rc <- sapply(1:nrow(chromsubs),function(i) flip(chromsubs$trinuc_Ref[i]))
        chromsubs$trinuc_Alt_rc <- sapply(1:nrow(chromsubs),function(i) flip(chromsubs$trinuc_Alt[i]))
        chromsubs$trinuc_mutation_rc <- paste0(chromsubs$trinuc_Ref_rc,">",chromsubs$trinuc_Alt_rc)
        # now check if the mutation is in the set of mutation types, otherwise use reverse complement
        isInMutTypes <- chromsubs$trinuc_Ref %in% chosen_trinucl
        chromsubs$trinuc_mutation_final <- sapply(1:nrow(chromsubs),function(x) {
          ifelse(isInMutTypes[x],chromsubs$trinuc_mutation[x],chromsubs$trinuc_mutation_rc[x])
        })
        sample_TNVs <- rbind(sample_TNVs,chromsubs)
      }
    }
    # add to final DNV table
    TNV_table <- rbind(TNV_table,sample_TNVs)
    
    #  now I have the sample DNVs, build also the catalogue
    if (is.null(sample_TNVs)){
      TNV_catalogue[,s] <- 0
    }else{
      countmuts <- table(sample_TNVs$trinuc_mutation_final)
      TNV_catalogue[,s] <- 0
      TNV_catalogue[names(countmuts),s] <- countmuts
    }
    
  }
  
  res <- list()
  res$TNV_catalogue <- TNV_catalogue
  res$TNV_table <- TNV_table
  return(res)
}

#' Build a Trinucleotide Variants Catalogue from TNVs list
#' 
#' This function takes as input a list of trinucleotide variants (TNVs). It then
#' annotates the TNVs and computes the TNV catalogues.
#' The trinucleotide mutations considered are only those where all three bases change,
#' i.e. TTT>AAA but not TTT>ATA.
#' 
#' @param tnvtab requires columns Sample, Chrom, Pos, Ref, Alt, with Ref and Alt of length 3
#' @return list of annotated TNVs and TNV catalogue
#' @export
tnvTabToTNVcatalogue <- function(tnvtab){
  
  # channels
  bases <- c("A","C","G","T")
  combos <- c()
  for (b1 in bases){
    for (b2 in bases){
      for (b3 in bases){
        combos <- c(combos,paste0(b1,b2,b3))
      }
    }
  }
  combos_rc <- sapply(combos,flip,USE.NAMES = F)
  # choose complement by alphabetical order
  chosen_trinucl <- unique(sapply(1:length(combos),function(i)ifelse(combos[i]>combos_rc[i],combos[i],combos_rc[i]),USE.NAMES = F))
  # now to get all mutation types
  mutationTypes <- c()
  for (trinucl in chosen_trinucl) {
    t1 <- substr(trinucl,1,1)
    t2 <- substr(trinucl,2,2)
    t3 <- substr(trinucl,3,3)
    for (possibleMut in combos) {
      # check if any base is the same
      m1 <- substr(possibleMut,1,1)
      m2 <- substr(possibleMut,2,2)
      m3 <- substr(possibleMut,3,3)
      if(t1 != m1 & t2 != m2 & t3 != m3) mutationTypes <- c(mutationTypes,paste(trinucl,possibleMut,sep = ">"))
    }
  }
  
  sample_list <- unique(tnvtab$Sample)
  TNV_catalogue <- data.frame(row.names = mutationTypes,stringsAsFactors = F)
  TNV_table <- NULL
  for (s in sample_list){
    # for debug: s <- sample_list[1]
    sample_TNVs <- NULL
    sample_TNVs <- tnvtab[tnvtab$Sample==s,,drop=F]
    sample_TNVs$trinuc_mutation <- paste0(sample_TNVs$Ref,">",sample_TNVs$Alt)
    # get reverse complement
    sample_TNVs$Ref_rc <- sapply(1:nrow(sample_TNVs),function(i) flip(sample_TNVs$Ref[i]))
    sample_TNVs$Alt_rc <- sapply(1:nrow(sample_TNVs),function(i) flip(sample_TNVs$Alt[i]))
    sample_TNVs$trinuc_mutation_rc <- paste0(sample_TNVs$Ref_rc,">",sample_TNVs$Alt_rc)
    # now check if the mutation is in the set of mutation types, otherwise use reverse complement
    isInMutTypes <- sample_TNVs$trinuc_mutation %in% mutationTypes
    sample_TNVs$trinuc_mutation_final <- sapply(1:nrow(sample_TNVs),function(x) {
      ifelse(isInMutTypes[x],sample_TNVs$trinuc_mutation[x],sample_TNVs$trinuc_mutation_rc[x])
    })
    
    # remove if not true TNV, for example if TAT>AAA, the central nucleotide is not changed
    sample_TNVs <- sample_TNVs[sample_TNVs$trinuc_mutation_final %in% mutationTypes,,drop=F]
    
    # add to final DNV table
    TNV_table <- rbind(TNV_table,sample_TNVs)
    
    #  now I have the sample DNVs, build also the catalogue
    if (is.null(sample_TNVs)){
      TNV_catalogue[,s] <- 0
    }else{
      countmuts <- table(sample_TNVs$trinuc_mutation_final)
      TNV_catalogue[,s] <- 0
      TNV_catalogue[names(countmuts),s] <- countmuts
    }
    
  }
  
  res <- list()
  res$TNV_catalogue <- TNV_catalogue
  res$TNV_table <- TNV_table
  return(res)
}

#' Plot Tinucleotide Variant Signatures or Catalogues
#' 
#' Function to plot one or more TNV signatures or catalogues. 
#' 
#' @param signature_data_matrix matrix of signatures or catalogues, signatures as columns and channels as rows.
#' @param output_file set output file, should end with ".jpg" of ".pdf". If output_file==null, output will not be to a file, but will still run the plot functions. The option output_file==null can be used to add this plot to a larger output file.
#' @param plot_sum whether the sum of the channels should be plotted. If plotting signatures this should be FALSE, but if plotting sample catalogues, this can be set to TRUE to display the number of mutations in each sample.
#' @param overall_title set the overall title of the plot
#' @param mar set the margin of the plot
#' @param howManyInOnePage how many signatures or catalogues should be plotted on one page. Multiple pages are plotted if more signatures/catalogues to plot have been requested
#' @param ncolumns how many columns should be used to arrange the signatures/catalogues to plot
#' @export
plotTNVcatalogues <- function(signature_data_matrix,
                              output_file = NULL,
                              plot_sum = TRUE,
                              overall_title = "",
                              add_to_titles = NULL,
                              mar = NULL,
                              howManyInOnePage=100,
                              ncolumns = 1){
  if(!is.null(output_file)) plottype <- substr(output_file,nchar(output_file)-2,nchar(output_file))
  # global parameters
  kelly_colors <- c('F3C300', '875692', 'F38400', 'A1CAF1', 'BE0032', 
                    'C2B280', '848482', '008856', 'E68FAC', '0067A5', 'F99379', '604E97', 
                    'F6A600', 'B3446C', 'DCD300', '882D17', '8DB600', '654522', 'E25822', '2B3D26','222222','F2F3F4', 'CCCCCC','CCCCCC','CCCCCC')
  kelly_colors <- paste0("#",kelly_colors)
  # repeatcolours <- c("red","green","orange","blue","purple","yellow")
  repeatcolours <- kelly_colors[1:6]
  
  signature_data_matrix <- as.matrix(signature_data_matrix[apply(signature_data_matrix,1,sum)>0,,drop=F])
  
  npages <- ceiling(ncol(signature_data_matrix)/howManyInOnePage)
  if(!is.null(output_file)) rootoutput_file <- substr(output_file,1,nchar(output_file)-4)
  for(i in 1:npages){
    ifrom <- howManyInOnePage*(i-1) + 1
    ito <- min(ncol(signature_data_matrix),howManyInOnePage*i)
    tmpmatrix <- signature_data_matrix[,ifrom:ito,drop=F]
    if (!is.null(add_to_titles)) tmpadd <- add_to_titles[ifrom:ito]
    if(npages>1 & !is.null(output_file)) output_file <- paste0(rootoutput_file,"_",i,"of",npages,".",plottype)
    nplotrows <- ceiling(ncol(tmpmatrix)/ncolumns)
    if(!is.null(output_file)) {
      if(plottype=="jpg"){
        jpeg(output_file,width = ncolumns*800,height = nplotrows*300,res = 220)
      }else if (plottype=="pdf"){
        pdf(output_file,width = ncolumns*8,height = nplotrows*3+0.5,pointsize = 26)
      }
      par(mfrow = c(nplotrows, ncolumns),omi=c(0,0,0.5,0),cex=0.7)
    }
    # keep only channels > 0, share same channels across samples
    if(nrow(tmpmatrix)>0){
      for (pos in 1:ncol(tmpmatrix)){
        if(is.null(mar)){
          par(mar=c(2,3.5,2,1))
        }else{
          par(mar=mar)
        }
        title <- colnames(tmpmatrix)[pos]
        if (plot_sum) title <- paste0(title," (",round(sum(tmpmatrix[,pos]))," TNVs)")
        if (!is.null(add_to_titles)) title <- paste0(title,"\n",tmpadd[pos])
        barplot(tmpmatrix[,pos],beside = T,border = NA,names.arg = "",col = repeatcolours,main = title)
      }
    }else{
      plot.new()
    }
    title(main = overall_title,outer = TRUE,cex.main = 1.5)
    if(!is.null(output_file)) dev.off()
  }
}

