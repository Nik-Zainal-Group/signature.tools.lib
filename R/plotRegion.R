#' plotRegion
#' 
#' Plotting of somatic variants in a specified region
#' 
#' @param outfilename file where the figure should be plotted. This also determines the file type of the output, use either pdf (recommended) or png
#' @param sample_name name of sample
#' @param chr chromosome of the region to plot
#' @param rstart start position of the region to plot
#' @param rend end position of the region to plot
#' @param plotSNVandIndels if TRUE then plot SNV and Indels at the bottom of the plot
#' @param SNV_vcf_file name of the vcf file containing the SNVs
#' @param SNV_tab_file name of the tab separated file containing the SNVs. Column names should be: chr, position, REF and ALT. If SNV_vcf_file is also specified, these variants will be ignored and the variants in the vcf file will be used instead.
#' @param Indels_vcf_file name of the vcf file containing the Indels
#' @param Indels_tab_file name of the tab separated file containing the Indels. Column names should be: chr, position, REF and ALT. If Indels_vcf_file is also specified, these variants will be ignored and the variants in the vcf file will be used instead.
#' @param CNV_tab_file name of the tab separated file containing the CNVs. Column names should be: 'seg_no', 'Chromosome', 'chromStart', 'chromEnd', 'total.copy.number.inNormal', 'minor.copy.number.inNormal', 'total.copy.number.inTumour', 'minor.copy.number.inTumour'
#' @param SV_bedpe_file name of the tab separated bedpe file containing the SVs. The file should contain a rearrangement for each row (two breakpoint positions should be on one row as determined by a pair of mates of paired-end sequencing) and should already be filtered according to the user preference, as all rearrangements in the file will be used and no filter will be applied. The files should contain a header in the first line with the following columns: "chrom1", "start1", "end1", "chrom2", "start2", "end2" and "sample" (sample name). In addition, either two columns indicating the strands of the mates, "strand1" (+ or -) and "strand2" (+ or -), or one column indicating the structural variant class, "svclass": translocation, inversion, deletion, tandem-duplication. The column "svclass" should correspond to (Sanger BRASS convention): inversion (strands +/- or -/+ and mates on the same chromosome), deletion (strands +/+ and mates on the same chromosome), tandem-duplication (strands -/- and mates on the same chromosome), translocation (mates are on different chromosomes). In addition, columns 'non-template'	and 'micro-homology' can be specified, including the non-templated insertion or micro-homology deletion sequence, which will be used to build an SV junctions catalogue.
#' @param plot_title title of the plot. If NULL, then the sample_name will be used as title. Use "", the empty string, to have no title.
#' @param genome.v genome version to use, either hg19 or hg38
#' @param spreadSVlabels when one of the two structural variant breakpoints is not in the region, it's location is written as text above the location of the breakpoint that is in the region. If spreadSVlabels is TRUE, all these text labels will be spread out evenly, which is useful if a lot of labels are close to each other and overlapping.
#' @param debug if TRUE, show debug guidelines and grid when plotting (default is FALSE)
#' @return all computed results, like catalogues and clustering regions, will be returned
#' @export
plotRegion <- function(outfilename,
                       sample_name,
                       chr,
                       rstart,
                       rend,
                       plotSNVandIndels = TRUE,
                       SNV_vcf_file = NULL,
                       SNV_tab_file = NULL,
                       snvs_table = NULL,
                       Indels_vcf_file = NULL,
                       Indels_tab_file = NULL,
                       indels_table = NULL,
                       CNV_tab_file = NULL,
                       CNV_table = NULL,
                       SV_bedpe_file = NULL,
                       sv_table = NULL,
                       plot_title = NULL,
                       spreadSVlabels = TRUE,
                       genome.v = "hg19",
                       debug=FALSE){
  
  # set up some variables
  # snvs_table <- NULL
  # indels_table <- NULL
  # CNV_table <- NULL
  # sv_table <- NULL
  
  sbs_obj <- NULL
  dbs_obj <- NULL
  indels_obj <- NULL
  sv_obj <- NULL
  
  # Loading SNVs if available
  if(!is.null(SNV_vcf_file)){
    if (file.exists(SNV_vcf_file)){
      if(is.null(snvs_table)){
        snvs_table <- fromVcfToTable(vcfFilename = SNV_vcf_file,
                                     genome.v = genome.v)
      }else{
        message("[warning plotRegion] SNV_vcf_file ignored because SNVs already loaded from snvs_table input parameter.")
        SNV_vcf_file <- NULL
      }
    }else{
      message("[warning plotRegion] SNV_vcf_file file not found: ",SNV_vcf_file,". Ignoring and moving on.")
      SNV_vcf_file <- NULL
    }
  }
  if(!is.null(SNV_tab_file)){
    if (file.exists(SNV_tab_file)){
      if(is.null(snvs_table)){
        snvs_table <- read.table(file = SNV_tab_file,sep = "\t",header = TRUE,
                                 check.names = FALSE,stringsAsFactors = FALSE)
      }else{
        message("[warning plotRegion] SNV_tab_file ignored because SNVs already loaded from SNV_vcf_file or from snvs_table input parameter.")
        SNV_tab_file <- NULL
      }
    }else{
      message("[warning plotRegion] SNV_tab_file file not found: ",SNV_tab_file,". Ignoring and moving on.")
      SNV_tab_file <- NULL
    }
  }
  if(!is.null(snvs_table)){
    # select only mutations in the given location
    lselection <- snvs_table$chr==chr & snvs_table$position >= rstart & snvs_table$position <= rend
    snvs_table <- snvs_table[lselection,,drop=F]
  }
  
  # Loading Indels if available
  if(!is.null(Indels_vcf_file)){
    if (file.exists(Indels_vcf_file)){
      # indels_obj <- vcfToIndelsClassification(indelsVCF.file = Indels_vcf_file,
      #                                         sampleID = sample_name,
      #                                         genome.v = genome.v)
      indels_table <- fromVcfToTable(vcfFilename = Indels_vcf_file,
                                     genome.v = genome.v)
    }else{
      message("[warning plotRegion] Indels_vcf_file file not found: ",Indels_vcf_file,". Ignoring and moving on.")
      Indels_vcf_file <- NULL
    }
  }
  if(!is.null(Indels_tab_file)){
    if (file.exists(Indels_tab_file)){
      if(is.null(indels_table)){
        indels_table <- read.table(file = Indels_tab_file,sep = "\t",header = TRUE,
                                   check.names = FALSE,stringsAsFactors = FALSE)
      }else{
        message("[warning plotRegion] Indels_tab_file ignored because Indels already loaded from Indels_vcf_file")
        Indels_tab_file <- NULL
      }
    }else{
      message("[warning plotRegion] Indels_tab_file file not found: ",Indels_tab_file,". Ignoring and moving on.")
      Indels_tab_file <- NULL
    }
  }
  if(!is.null(indels_table)){
    # select only mutations in the given location
    lselection <- indels_table$chr==chr & indels_table$position >= rstart & indels_table$position <= rend
    indels_table <- indels_table[lselection,,drop=F]
  }
  
  # Loading CNV file
  if(!is.null(CNV_tab_file)){
    if (file.exists(CNV_tab_file)){
      CNV_table <- read.table(file = CNV_tab_file,sep = "\t",header = TRUE,
                              check.names = FALSE,stringsAsFactors = FALSE)
    }else{
      message("[warning plotRegion] CNV_tab_file file not found: ",CNV_tab_file,". Ignoring and moving on.")
      CNV_tab_file <- NULL
    }
  }
  if(!is.null(CNV_table)){
    CNV_table <- slice_cnv_table_to_regions(cnv_table = CNV_table,
                                            regionsDF = data.frame(chr=chr,
                                                                   start=rstart,
                                                                   end=rend,
                                                                   stringsAsFactors = F))
  }
  
  # Loading SV bedpe file
  if(!is.null(SV_bedpe_file)){
    if (file.exists(SV_bedpe_file)){
      sv_table <- readTable(file = SV_bedpe_file)
      
      #check whether column svclass is present,
      #if not, compute it
      if (! "svclass" %in% colnames(sv_table)){
        if ("strand1" %in% colnames(sv_table) & "strand2" %in% colnames(sv_table)){
          sv_table <- classifyRearrangementsFromBedpe(sv_table)
        }else{
          message("[warning plotRegion] cannot classify structural variants: svclass column missing, and cannot compute it because strand1 and strand2 are missing.")
        }
      }
      
    }else{
      message("[warning plotRegion] SV_bedpe_file file not found: ",SV_bedpe_file,". Ignoring and moving on.")
      SV_bedpe_file <- NULL
    }
  }
  if(!is.null(sv_table)){
    sv_table <- slice_sv_table_to_regions(sv_table = sv_table,
                                          regionsDF = data.frame(chr=chr,
                                                                 start=rstart,
                                                                 end=rend,
                                                                 stringsAsFactors = F))
  }
  
  # check if there is anything at all to plot
  if(is.null(snvs_table) & is.null(indels_table) & is.null(CNV_table) & is.null(sv_table)){
    message("[error plotRegion] no data found, nothing to plot. Quit.")
    return(NULL)
  }
  
  # get SBS catalogue if possible
  snvs_classified <- NULL
  if(!is.null(snvs_table)){
    if(nrow(snvs_table)>0){
      sbs_obj <- tabToSNVcatalogue(subs = snvs_table,
                                   genome.v = genome.v)
      colnames(sbs_obj$muts)[c(1,2,4,5)] <- c("chr","position","REF","ALT")
      snvs_classified <- calcIntermutDist(sbs_obj$muts)
      colnames(sbs_obj$catalogue) <- "SNV catalogue"
    }
  }
  
  # get indels classification if possible
  if(!is.null(indels_table)){
    if(nrow(indels_table)>0){
      indels_obj <- tabToIndelsClassification(indel.data = indels_table,
                                              sampleID = sample_name,
                                              genome.v = genome.v)
    }
  }
  
  # get DBS catalogue if possible
  if(!is.null(snvs_table)){
    if(nrow(snvs_table)>1){
      dbs_obj <- tabToDNVcatalogue(muttable = snvs_table)
      colnames(dbs_obj$DNV_catalogue) <- "DNV catalogue"
    }
  }
  
  # get sv classification if possible
  if(!is.null(sv_table)){
    if(nrow(sv_table)>0){
      sv_obj <- bedpeToRearrCatalogue(sv_bedpe = sv_table)
      colnames(sv_obj$rearr_catalogue) <- "SV catalogue"
      if(!is.null(sv_obj$junctions_catalogue)) colnames(sv_obj$junctions_catalogue) <- "SV junction catalogue"
    }
  }
  
  # time to plot, outfilename needs to be specified
  plottype <- substr(outfilename, nchar(outfilename) - 2, nchar(outfilename))
  dir.create(dirname(outfilename),showWarnings = F,recursive = T)
  
  # function for debug of plots position
  drawDebugBox <- function(n){
    kelly_colors <- getKellyColours()
    
    par(mai=c(0,0,0,0),new=TRUE)
    plot(NULL,xlim = c(0,1),ylim = c(0,1),main="",type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",xaxs="i",yaxs="i")
    rect(0,0,1,1,border = kelly_colors[n],lwd = 3,lty = 3)
  }
  
  # organise panels according to a certain plot dimension (x,y)
  plotx <- 14
  if(plotSNVandIndels){
    ploty <- 6.5
  }else{
    ploty <- 6.5 - 1.75
  }
  
  debug_grid_gap <- 1
  debug_grid_tick <- 0.05
  msgtextscaling <- 0.5
  
  placePanel <- function(where,width,height){
    par(fig = c(where[1]/plotx,
                (where[1]+width)/plotx,
                where[2]/ploty,
                (where[2]+height)/ploty),
        new=TRUE)
  }
  
  # open the file
  if(plottype=="pdf"){
    cairo_pdf(filename = outfilename,width = plotx,height = ploty)
  }else if(plottype=="png"){
    png(filename = outfilename,width = 300*plotx,height = 300*ploty,res = 300)
  }else{
    message("[error plotRegion] incorrect file type: ",plottype,". ",
            "Please end your file name with .pdf or .png")
    return(NULL)
  }
  
  # plot data
  if(is.null(plot_title)) plot_title <- sample_name
  if(debug){
    par(fig = c(0,1,0,1),mai=c(0,0,0,0))
    plot(NULL,xlim = c(0,plotx),ylim = c(0,ploty),main="",type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",xaxs="i",yaxs="i")
    gridx <- seq(0,plotx,debug_grid_gap)
    gridy <- seq(0,ploty,debug_grid_gap)
    abline(v = gridx,col="grey",lty = 3)
    abline(h = gridy,col="grey",lty = 3)
    tickx <- seq(0,plotx,debug_grid_tick)
    ticky <- seq(0,ploty,debug_grid_tick)
    for(gy in gridy) for(xi in tickx) lines(c(xi,xi),c(gy-0.004,gy+0.004),col="grey")
    for(gx in gridx) for(yi in ticky) lines(c(gx-0.004,gx+0.004),c(yi,yi),col="grey")
    text(x = gridx[1:(length(gridx)-1)]+0.007,y = 0.01,
         labels = gridx[1:(length(gridx)-1)],col = "grey",cex=0.5,adj = 0)
    text(x = 0.007,y = gridy[1:(length(gridy)-1)]+0.01,
         labels = gridy[1:(length(gridy)-1)],col = "grey",cex=0.5,adj = 0)
    par(fig = c(0,1,0,1),mai=c(0,0.1,0.4,0),new=TRUE)
  }else{
    par(fig = c(0,1,0,1),mai=c(0,0.1,0.4,0))
  }
  plot(NULL,xlim = c(0,1),ylim = c(0,1),main="",type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",xaxs="i",yaxs="i")
  title(main = plot_title,adj = 0)
  
  # SNVs colours
  snvcolours <- getSBSCatalogueColours()
  # chromosome label
  xlab <- ifelse(startsWith(as.character(chr),prefix = "chr"),substr(chr,4,8),chr)
  
  # plot the SNVs at the bottom
  if(plotSNVandIndels){
    placePanel(where = c(0,0),width = 14,height =6.5)
    par(mai=c(1,1,4,0.5))
    plot(x = snvs_classified$position[-1]/1e6,
         y = log10(snvs_classified$distPrev[-1]),
         xlim = c(rstart/1e6,rend/1e6),
         ylim = c(0,8),
         las=1,
         xaxs="i",
         pch=16,
         col=snvcolours[paste0(snvs_classified$pyrwt,">",snvs_classified$pyrmut)],
         xlab = paste0("chromosome ",xlab," (Mb)"),
         ylab = "log10(IMD)")
    if(debug) drawDebugBox(1)
    
    #plot indels
    indels_colours <- getIndelsClassColours()
    placePanel(where = c(0,2.5),width = 14,height =0.25)
    par(mai=c(0,1,0,0.5))
    plot(x = NULL,
         xlim = c(rstart,rend),
         ylim = c(0,1),
         xaxs="i",
         xaxt="n",
         yaxt="n",
         xlab = "",
         ylab = "")
    # abline(h=1)
    mtext(text = "indels",side=2,las=1,line=2.2)
    if(!is.null(indels_obj)){
      if(nrow(indels_obj$indels_classified)>0){
        for (i in 1:nrow(indels_obj$indels_classified)) {
          # i <- 1
          lines(x=c(indels_obj$indels_classified$pos[i],
                    indels_obj$indels_classified$pos[i]),
                y=c(0,1),
                lwd=2,
                col=indels_colours[indels_obj$indels_classified$indel.class[i]])
        }
      }
    }
    if(debug) drawDebugBox(1)
  }

  
  # plot copy numbers
  max_CN <- 2
  
  if(!is.null(CNV_table)){
    if(nrow(CNV_table)>0){
      max_CN <- max(max_CN,CNV_table$total.copy.number.inTumour)
    }
  }
  
  totalCNcolour <- "#9ACD32E6"
  
  if(plotSNVandIndels){
    placePanel(where = c(0,2.75),width = 14,height = 1)
    par(mai=c(0,1,0,0.5))
    plot(x = NULL,
         xlim = c(rstart,rend)/1e6,
         ylim = c(0,max_CN),
         las=1,
         xaxs="i",
         xaxt="n",
         xlab = "",
         ylab = "Total CN")
  }else{
    placePanel(where = c(0,0),width = 14,height = 2)
    par(mai=c(1,1,0,0.5))
    plot(x = NULL,
         xlim = c(rstart,rend)/1e6,
         ylim = c(0,max_CN),
         las=1,
         xaxs="i",
         xlab = paste0("chromosome ",xlab," (Mb)"),
         ylab = "Total CN")
  }

  
  if(!is.null(CNV_table)){
    if(nrow(CNV_table)>0){

      for (i in 1:nrow(CNV_table)) {
        # i <- 1
        lines(x = c(CNV_table$chromStart[i],
                    CNV_table$chromEnd[i])/1e6,
              y=c(CNV_table$total.copy.number.inTumour[i],
                  CNV_table$total.copy.number.inTumour[i]),
              col=totalCNcolour,
              lwd=2)
      }
    }
  }
  
  # plot structural variants 
  # placePanel(where = c(0,0),width = 14,height =7)
  # par(mai=c(1,1,4,0.5))
  if(plotSNVandIndels){
    placePanel(where = c(0,2.75),width = 14,height = 3.75)
  }else{
    placePanel(where = c(0,1),width = 14,height = 3.75)
  }
  par(mai=c(0,1,2.75,0.5))
  plot(x = NULL,
       xlim = c(rstart,rend),
       ylim = c(0,4),
       xaxs="i",
       type="n",
       xaxt="n",
       yaxt="n",
       xlab="",
       ylab="",
       bty="n")
  if(debug) drawDebugBox(2)
  
  # plotting structural variants if any
  # moveup <- 5.5
  moveup <- 4
  # SVcol <- getSVClassColours()
  SVcol <- getSVClassColours_SVCatalogue()
  par(xpd=T)
  # ysvpos <- 8+moveup+c(1:3,5)
  ysvpos <- moveup+c(1:3,5)
  names(ysvpos) <- c("tandem-duplication",
                     "deletion",
                     "inversion",
                     "translocation")
  for(p in ysvpos){
    lines(x=c(rstart,rend),
          y=c(p,p),
          lty=2,
          lwd=2,
          col="grey")
  }
  
  # get a new order for plotting, so that I can plot the SVs that have only one
  # breakpoint in the region first, and also I can plot them in order and assign
  # an x position to each label for the other breakpoint
  singletonsRows <- c()
  singletonsXpos <- c()
  if(nrow(sv_table)>0){
    for(svi in 1:nrow(sv_table)){
      # svi <- 1
      if(sv_table$bp1inregion[svi] & !sv_table$bp2inregion[svi]){
        # get bp1 pos
        xpos <- (sv_table$start1[svi]+sv_table$end1[svi])/2
        singletonsRows <- c(singletonsRows,svi)
        singletonsXpos <- c(singletonsXpos,xpos)
      }else if(!sv_table$bp1inregion[svi] & sv_table$bp2inregion[svi]){
        # get bp2 pos
        xpos <- (sv_table$start2[svi]+sv_table$end2[svi])/2
        singletonsRows <- c(singletonsRows,svi)
        singletonsXpos <- c(singletonsXpos,xpos)
      }
    }
  }
  
  if(length(singletonsRows)>0){
    newroworder <- c(singletonsRows[order(singletonsXpos)],setdiff(1:nrow(sv_table),singletonsRows))
    xlabelpos <- singletonsXpos[order(singletonsXpos)]
    xlabeli <- 1
    if(spreadSVlabels){
      xlabelpos <- seq(rstart,rend,length.out=length(singletonsRows))
    }
  }else{
    if(nrow(sv_table)>0) newroworder <- 1:nrow(sv_table)
  }

  
  if(nrow(sv_table)>0){
    for(svi in newroworder){
      # svi <- 1
      if(sv_table$bp1inregion[svi] & sv_table$bp2inregion[svi]){
        # need to plot them together and connect them if possible
        ypos <- ysvpos[sv_table$svclass[svi]]
        xpos1 <- (sv_table$start1[svi]+sv_table$end1[svi])/2
        lines(x=c(xpos1,xpos1),
              y=c(0,ypos),
              lwd=2,
              col=SVcol[sv_table$svclass[svi]])
        xpos2 <- (sv_table$start2[svi]+sv_table$end2[svi])/2
        lines(x=c(xpos2,xpos2),
              y=c(0,ypos),
              lwd=2,
              col=SVcol[sv_table$svclass[svi]])
        # draw arc
        n <- 100
        b <- 1
        a <- (xpos2-xpos1)/2
        xarc <- seq(from = 0,to = a, length.out = n)
        yarc <- sapply(xarc,function(x){
          b/a*sqrt(a*a-x*x)
        })
        lines(x=xpos1+a+xarc,
              y=ypos+yarc,
              lwd=2,
              col=SVcol[sv_table$svclass[svi]])
        lines(x=xpos1+a-xarc,
              y=ypos+yarc,
              lwd=2,
              col=SVcol[sv_table$svclass[svi]])
      }else if(sv_table$bp1inregion[svi]){
        # plot bp1
        xpos <- (sv_table$start1[svi]+sv_table$end1[svi])/2
        lines(x=c(xpos,xpos),
              y=c(0,ysvpos[sv_table$svclass[svi]]),
              lwd=2,
              col=SVcol[sv_table$svclass[svi]])
        # where is the other bp?
        xchr2 <- sv_table$chrom2[svi]
        xpos2 <- (sv_table$start2[svi]+sv_table$end2[svi])/2
        xlp <- xlabelpos[xlabeli]
        xlabeli <- xlabeli+1
        
        lines(x=c(xpos,xlp),
              y=c(max(ysvpos)+0.2,max(ysvpos)+1))
        text(x=xlp,
             y=max(ysvpos)+1.2,
             labels = paste0(xchr2,":",format(floor(xpos2), big.mark=",", scientific=FALSE)),
             srt=90,
             adj = 0,
             cex = 0.6)
      }else if(sv_table$bp2inregion[svi]){
        # plot bp2
        xpos <- (sv_table$start2[svi]+sv_table$end2[svi])/2
        lines(x=c(xpos,xpos),
              y=c(0,ysvpos[sv_table$svclass[svi]]),
              lwd=2,
              col=SVcol[sv_table$svclass[svi]])
        # where is the other bp?
        xchr2 <- sv_table$chrom1[svi]
        xpos2 <- (sv_table$start1[svi]+sv_table$end1[svi])/2
        xlp <- xlabelpos[xlabeli]
        xlabeli <- xlabeli+1
        
        lines(x=c(xpos,xlp),
              y=c(max(ysvpos)+0.2,max(ysvpos)+1))
        text(x=xlp,
             y=max(ysvpos)+1.2,
             labels = paste0(xchr2,":",format(floor(xpos2), big.mark=",", scientific=FALSE)),
             srt=90,
             adj = 0,
             cex = 0.6)
      }
    }
  }
  par(xpd=F)

  
  # close the file
  dev.off()
  
  # also I can return all the data/info obtained
  returnObj <- NULL
  returnObj$sample_name <- sample_name
  returnObj$genomev <- genome.v
  returnObj$snvs_table <- snvs_table
  returnObj$indels_table <- indels_table
  returnObj$CNV_table <- CNV_table
  returnObj$sv_table <- sv_table
  returnObj$sbs_obj <- sbs_obj
  returnObj$dbs_obj <- dbs_obj
  returnObj$indels_obj <- indels_obj
  returnObj$sv_obj <- sv_obj
  # returnObj$kataegis_regions <- kataegis_regions
  # returnObj$kataegisSBScatalogue_all <- kataegisSBScatalogue_all
  # returnObj$clustering_regions_sbs_catalogues <- clustering_regions_sbs_catalogues
  # returnObj$clusteringSBScatalogue_all <- clusteringSBScatalogue_all
  return(returnObj)
  
}

slice_cnv_table_to_regions <- function(cnv_table,
                                       regionsDF){
  tmp_table <- NULL
  for(i in 1:nrow(regionsDF)){
    # i <- 1
    selection <- cnv_table$Chromosome==regionsDF$chr[i]
    tmprows <- cnv_table[selection,,drop=F]
    if(nrow(tmprows)>0){
      selection <- !(tmprows$chromStart > regionsDF$end[i] | tmprows$chromEnd < regionsDF$start[i])
      tmprows <- tmprows[selection,,drop=F]
      if(nrow(tmprows)>0){
        # slice off
        tmprows$chromStart[tmprows$chromStart<regionsDF$start[i]] <- regionsDF$start[i]
        tmprows$chromEnd[tmprows$chromEnd>regionsDF$end[i]] <- regionsDF$end[i]
        tmp_table <- rbind(tmp_table,tmprows)
      }
    }
  }
  return(tmp_table)
}

# bothBPinRegion: if TRUE we keep only the sv with both BP in the region, 
#                 otherwise only one BP is enough
slice_sv_table_to_regions <- function(sv_table,
                                      regionsDF,
                                      bothBPinRegion=F){
  selection1 <- rep(FALSE,nrow(sv_table))
  selection2 <- rep(FALSE,nrow(sv_table))
  for(i in 1:nrow(regionsDF)){
    # i <- 1
    selection1 <- selection1 | (sv_table$chrom1==regionsDF$chr[i] & !(sv_table$start1 > regionsDF$end[i] | sv_table$end1 < regionsDF$start[i]))
    selection2 <- selection2 | (sv_table$chrom2==regionsDF$chr[i] & !(sv_table$start2 > regionsDF$end[i] | sv_table$end2 < regionsDF$start[i]))
  }
  # add the annotation
  sv_table$bp1inregion <- selection1
  sv_table$bp2inregion <- selection2
  if(bothBPinRegion){
    selection <- selection1 & selection2
  }else{
    selection <- selection1 | selection2
  }
  return(sv_table[selection,,drop=F])
}

