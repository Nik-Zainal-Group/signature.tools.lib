#' Transcription and Replication strand bias for a sample
#'
#' This function computes the transcription/replication strand bias of a given
#' set of single base substitutions, typically from a single sample. The function
#' returns the bias of the mutations with and without the trinucleotide context.
#'
#' @param snv_table data frame with subs from a single sample and the following minimal columns: chr, position, REF, ALT.
#' @param genome.v either "hg38" (will load BSgenome.Hsapiens.NCBI.GRCh38), "hg19" (will load BSgenome.Hsapiens.UCSC.hg19), mm10 (will load BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10) or canFam3 (will load BSgenome.Cfamiliaris.UCSC.canFam3::BSgenome.Cfamiliaris.UCSC.canFam3)
#' @return strand bias results
#' @export
sampleStrandBias <- function(snv_table,
                             genomev = "hg19"){
  
  # checking for necessary input columns
  if(!all(c("chr","position","REF","ALT") %in% colnames(snv_table))){
    message("[error sampleStrandBias] missing columns in the input snv_table.",
            "Please make sure the following column names are present: chr, position, REF, ALT.")
    return(NULL)
  }
  
  # set up some variables
  bias_types <- c("transcription","replication")
  bias_subtypes <- list()
  bias_subtypes[["transcription"]] <- c("uts","ts")
  bias_subtypes[["replication"]] <- c("leading","lagging")
  # some more variables
  all_bp <- c("A", "C", "G", "T")
  pyr_muts <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  pyr_trimuts <- c()
  for(m in pyr_muts){
    for(a in all_bp){
      for(b in all_bp){
        pyr_trimuts <- c(pyr_trimuts, paste(a, "[",m, "]", b, sep=""))
      }
    }
  }
  
  # set the genome version
  if(genomev=="hg19"){
    expected_chroms <- paste0(c(seq(1:22),"X","Y"))
    genomeSeq <- BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
  }else if(genomev=="hg38"){
    expected_chroms <- paste0("chr",c(seq(1:22),"X","Y"))
    genomeSeq <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  }
  
  # fetch is the expected bias
  strandBiasCounts_single_gv <- strandBiasCounts_single[strandBiasCounts_single$genomeversion==genomev,2:ncol(strandBiasCounts_single)]
  strandBiasCounts_tri_gv <- strandBiasCounts_tri[strandBiasCounts_tri$genomeversion==genomev,2:ncol(strandBiasCounts_tri)]
  
  # prepare output
  tmpTable <- unique(strandBiasTable[,2:3])
  rownames(tmpTable) <- 1:nrow(tmpTable)
  bias_results_single <- cbind(tmpTable,data.frame(matrix(0,
                                                          nrow = nrow(tmpTable),
                                                          ncol = length(pyr_muts),
                                                          dimnames = list(1:nrow(tmpTable),pyr_muts)),
                                                   stringsAsFactors = F,check.names = F))
  bias_results_tri <- cbind(tmpTable,data.frame(matrix(0,
                                                       nrow = nrow(tmpTable),
                                                       ncol = length(pyr_trimuts),
                                                       dimnames = list(1:nrow(tmpTable),pyr_trimuts)),
                                                stringsAsFactors = F,check.names = F))
  
  for (bt in bias_types) {
    # bt <- bias_types[1]
    for(bst in bias_subtypes[[bt]]) {
      # bst <- bias_subtypes[[bt]][1]
      tableSelection <- strandBiasTable$genomeversion==genomev & strandBiasTable$biastype==bt & strandBiasTable$biassubtype==bst
      table_slice <- strandBiasTable[tableSelection,,drop=F]
      # fix some stuff
      if(genomev=="hg19") table_slice$chr <- substr(table_slice$chr,4,10)
      table_slice <- table_slice[table_slice$chr %in% expected_chroms,,drop=F]
      # let's see which mutations are in these locations
      tmp_snv_table <- snv_table[sapply(1:nrow(snv_table),function(x) any(snv_table[x,"chr"]==table_slice$chr & snv_table[x,"position"]>=table_slice$start & snv_table[x,"position"]<=table_slice$end)),,drop=F]
      
      muts_single <- table(paste0(tmp_snv_table$REF,">",tmp_snv_table$ALT))
      muts_5p <- as.character(BSgenome::getSeq(genomeSeq,
                                               names=tmp_snv_table$chr,
                                               start=tmp_snv_table$position-1,
                                               end=tmp_snv_table$position-1))
      muts_3p <- as.character(BSgenome::getSeq(genomeSeq,
                                               names=tmp_snv_table$chr,
                                               start=tmp_snv_table$position+1,
                                               end=tmp_snv_table$position+1))
      muts_tri <- table(paste0(muts_5p,"[",tmp_snv_table$REF,">",tmp_snv_table$ALT,"]",muts_3p))
      
      # add the counts accordingly
      current_muts <- intersect(names(muts_single),pyr_muts)
      if(length(current_muts)>0) bias_results_single[bias_results_single$biastype==bt & bias_results_single$biassubtype==bst,current_muts] <- bias_results_single[bias_results_single$biastype==bt & bias_results_single$biassubtype==bst,current_muts] + muts_single[current_muts]
      current_muts <- setdiff(names(muts_single),pyr_muts)
      if(length(current_muts)>0) bias_results_single[bias_results_single$biastype==bt & bias_results_single$biassubtype==setdiff(bias_subtypes[[bt]],bst),sapply(current_muts,toPyrMutation,simplify = T,USE.NAMES = F)] <- bias_results_single[bias_results_single$biastype==bt & bias_results_single$biassubtype==setdiff(bias_subtypes[[bt]],bst),sapply(current_muts,toPyrMutation,simplify = T,USE.NAMES = F)] + muts_single[current_muts]
      
      current_muts <- intersect(names(muts_tri),pyr_trimuts)
      if(length(current_muts)>0) bias_results_tri[bias_results_single$biastype==bt & bias_results_single$biassubtype==bst,current_muts] <- bias_results_tri[bias_results_single$biastype==bt & bias_results_single$biassubtype==bst,current_muts] + muts_tri[current_muts]
      current_muts <- setdiff(names(muts_tri),pyr_trimuts)
      if(length(current_muts)>0) bias_results_tri[bias_results_single$biastype==bt & bias_results_single$biassubtype==setdiff(bias_subtypes[[bt]],bst),sapply(current_muts,toPyrMutation,simplify = T,USE.NAMES = F)] <- bias_results_tri[bias_results_single$biastype==bt & bias_results_single$biassubtype==setdiff(bias_subtypes[[bt]],bst),sapply(current_muts,toPyrMutation,simplify = T,USE.NAMES = F)] + muts_tri[current_muts]
      
    }
  }
  
  # return results
  returnObj <- list()
  returnObj$bias_results_single <- bias_results_single
  returnObj$bias_results_tri <- bias_results_tri
  returnObj$bias_results_single_ratios <- biasCountsToRatios(bias_results_single)
  returnObj$bias_results_tri_ratios <- biasCountsToRatios(bias_results_tri)
  returnObj$bias_expected_single_ratios <- biasCountsToRatios(strandBiasCounts_single_gv)
  returnObj$bias_expected_tri_ratios <- biasCountsToRatios(strandBiasCounts_tri_gv)
  
  return(returnObj)
}

toPyrMutation <- function(x){
  pyr_bases <- c("C","T")
  names(pyr_bases) <- c("G","A")
  
  bp <- c("A", "C", "G", "T")
  names(bp) <- c("T", "G", "C", "A")
  
  pyr_muts <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  names(pyr_muts) <- c("G>T", "G>C", "G>A", "A>T", "A>G", "A>C")
  
  if(nchar(x)==3){
    if(!(substr(x,1,1) %in% pyr_bases)){
      return(paste0(pyr_muts[x]))
    }else{
      return(x)
    }
  }else if(nchar(x)==7){
    if(!(substr(x,3,3) %in% pyr_bases)){
      return(paste0(bp[substr(x,7,7)],"[",pyr_muts[substr(x,3,5)],"]",bp[substr(x,1,1)]))
    }else{
      return(x)
    }
  }else{
    return(NULL)
  }
}

mutToBases <- function(x){
  if(nchar(x)==3){
    return(substr(x,1,1))
  }else if(nchar(x)==7){
    return(paste0(substr(x,1,1),substr(x,3,3),substr(x,7,7)))
  }else{
    return(NULL)
  }
}

biasCountsToRatios <- function(bias_counts_table){
  bias_types <- unique(bias_counts_table$biastype)
  features <- colnames(bias_counts_table)[3:ncol(bias_counts_table)]
  resultTable <- data.frame(features,
                            stringsAsFactors = F,
                            check.names = F)
  for(bt in bias_types){
    # bt <- bias_types[1]
    bias_subtypes <- unique(bias_counts_table$biassubtype[bias_counts_table$biastype==bt])
    resultTable[,paste0(bias_subtypes[1],"/",bias_subtypes[2])] <- unlist(bias_counts_table[bias_counts_table$biastype==bt & bias_counts_table$biassubtype==bias_subtypes[1],3:ncol(bias_counts_table)]/bias_counts_table[bias_counts_table$biastype==bt & bias_counts_table$biassubtype==bias_subtypes[2],3:ncol(bias_counts_table)])
  }
  return(resultTable)
}

#' plot Transcription and Replication strand bias results
#'
#' This function plots the transcription/replication strand bias results obtained
#' using the sampleStrandBias function.
#'
#' @param biasResObj R list object containing the results from the sampleStrandBias function call
#' @param filename pdf file name to save the plots to file
#' @param addToTitle text to be added in each plot title, useful for example to add the name of a sample
#' @export
plotBiasResults <- function(biasResObj,
                            filename=NULL,
                            addToTitle=NULL){
  # plot
  kelly_colors <- c('F3C300', '875692', 'F38400', 'A1CAF1', 'BE0032', 
                    'C2B280', '848482', '008856', 'E68FAC', '0067A5', 'F99379', '604E97', 
                    'F6A600', 'B3446C', 'DCD300', '882D17', '8DB600', '654522', 'E25822', '2B3D26','222222','F2F3F4', 'CCCCCC','CCCCCC','CCCCCC')
  kelly_colors <- paste0("#",kelly_colors)
  
  bias_ratios <- c("uts/ts","leading/lagging")
  biasTitles <- list()
  biasTitles[[bias_ratios[1]]] <- "Transcription bias"
  biasTitles[[bias_ratios[2]]] <- "Replication bias"
  gap <- max(abs(biasResObj$bias_results_single_ratios[,2:3] - 1))*1.2
  
  if(!is.null(filename)) cairo_pdf(filename = filename,width = 10,height = 5,pointsize = 12)
  par(mfrow=c(1,2),cex=1)
  par(mar=c(5,4,3,1))
  for(br in bias_ratios){
    usetitle <- biasTitles[[br]]
    if (!is.null(addToTitle)) usetitle <- paste0(usetitle,addToTitle)
    plot(x = biasResObj$bias_results_single_ratios[,br],
         y = nrow(biasResObj$bias_results_single_ratios):1,
         pch = 16,
         col = kelly_colors[3],
         yaxt='n',
         main= usetitle,
         xlim = c(1-gap,1+gap),xlab = paste0("ratio ",br),ylab = "")
    points(x=rep(biasResObj$bias_expected_single_ratios[,br],each=3),
           y=nrow(biasResObj$bias_results_single_ratios):1,
           col = kelly_colors[4],
           pch = 17)
    axis(side = 2,at = nrow(biasResObj$bias_results_single_ratios):1,
         labels = biasResObj$bias_results_single_ratios$features,las=1)
    abline(v=1,lty=2)
    legend(x="topright",legend = c("observed","expected"),fill = kelly_colors[3:4],border = F,bty = 'n')
  }
  if(!is.null(filename)) dev.off()
  
}
