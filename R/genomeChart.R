#' genomeChart
#' 
#' Plotting of somatic variants using the circlize R package. Somatic variants
#' plotted are single nucleotide variants (SNVs), Indels, copy number variants
#' (CNVs) and structural variants (SVs or rearrangements). In addition to the
#' circle visualisation of the variants, other data are plotted, such as SNV and
#' SV mutational catalogues, indels classification summary counts and a copy number
#' plot indicating estimates of total and minor allele copy numbers.
#' 
#' @param outfilename file where the figure should be plotted. This also determines the file type of the output, use either pdf (recommended) or png
#' @param sample_name name of sample
#' @param SNV_vcf_file name of the vcf file containing the SNVs
#' @param SNV_tab_file name of the tab separated file containing the SNVs. Column names should be: chr, position, REF and ALT. If SNV_vcf_file is also specified, these variants will be ignored and the variants in the vcf file will be used instead.
#' @param Indels_vcf_file name of the vcf file containing the Indels
#' @param Indels_tab_file name of the tab separated file containing the Indels. Column names should be: chr, position, REF and ALT. If Indels_vcf_file is also specified, these variants will be ignored and the variants in the vcf file will be used instead.
#' @param CNV_tab_file name of the tab separated file containing the CNVs. Column names should be: 'seg_no', 'Chromosome', 'chromStart', 'chromEnd', 'total.copy.number.inNormal', 'minor.copy.number.inNormal', 'total.copy.number.inTumour', 'minor.copy.number.inTumour'
#' @param SV_bedpe_file name of the tab separated bedpe file containing the SVs. The file should contain a rearrangement for each row (two breakpoint positions should be on one row as determined by a pair of mates of paired-end sequencing) and should already be filtered according to the user preference, as all rearrangements in the file will be used and no filter will be applied. The files should contain a header in the first line with the following columns: "chrom1", "start1", "end1", "chrom2", "start2", "end2" and "sample" (sample name). In addition, either two columns indicating the strands of the mates, "strand1" (+ or -) and "strand2" (+ or -), or one column indicating the structural variant class, "svclass": translocation, inversion, deletion, tandem-duplication. The column "svclass" should correspond to (Sanger BRASS convention): inversion (strands +/- or -/+ and mates on the same chromosome), deletion (strands +/+ and mates on the same chromosome), tandem-duplication (strands -/- and mates on the same chromosome), translocation (mates are on different chromosomes). In addition, columns 'non-template'	and 'micro-homology' can be specified, including the non-templated insertion or micro-homology deletion sequence, which will be used to build an SV junctions catalogue.
#' @param plot_title title of the plot. If NULL, then the sample_name will be used as title. Use "", the empty string, to have no title.
#' @param runKataegis whether or not to run the kataegis algorithm (default is TRUE)
#' @param genome.v genome version to use, either hg19 or hg38
#' @param debug if TRUE, show debug guidelines and grid when plotting (default is FALSE)
#' @return all computed results, like catalogues and clustering regions, will be returned
#' @export
genomeChart <- function(outfilename,
                        sample_name,
                        SNV_vcf_file = NULL,
                        SNV_tab_file = NULL,
                        Indels_vcf_file = NULL,
                        Indels_tab_file = NULL,
                        CNV_tab_file = NULL,
                        SV_bedpe_file = NULL,
                        plot_title = NULL,
                        runKataegis = TRUE,
                        genome.v = "hg19",
                        debug=FALSE){
  
  # set up some variables
  snvs_table <- NULL
  sbs_obj <- NULL
  dbs_obj <- NULL
  indels_obj <- NULL
  CNV_table <- NULL
  sv_obj <- NULL
  
  # Loading SNVs if available
  if(!is.null(SNV_vcf_file)){
    if (file.exists(SNV_vcf_file)){
      snvs_table <- fromVcfToTable(vcfFilename = SNV_vcf_file,
                                   genome.v = genome.v)
    }else{
      message("[warning genomeChart] SNV_vcf_file file not found: ",SNV_vcf_file,". Ignoring and moving on.")
      SNV_vcf_file <- NULL
    }
  }
  if(!is.null(SNV_tab_file)){
    if (file.exists(SNV_tab_file)){
      if(is.null(snvs_table)){
        snvs_table <- read.table(file = SNV_tab_file,sep = "\t",header = TRUE,
                                 check.names = FALSE,stringsAsFactors = FALSE)
      }else{
        message("[warning genomeChart] SNV_tab_file ignored because SNVs already loaded from SNV_vcf_file.")
        SNV_tab_file <- NULL
      }
    }else{
      message("[warning genomeChart] SNV_tab_file file not found: ",SNV_tab_file,". Ignoring and moving on.")
      SNV_tab_file <- NULL
    }
  }
  
  # Loading Indels if available
  if(!is.null(Indels_vcf_file)){
    if (file.exists(Indels_vcf_file)){
      indels_obj <- vcfToIndelsClassification(indelsVCF.file = Indels_vcf_file,
                                              sampleID = sample_name,
                                              genome.v = genome.v)
    }else{
      message("[warning genomeChart] Indels_vcf_file file not found: ",Indels_vcf_file,". Ignoring and moving on.")
      Indels_vcf_file <- NULL
    }
  }
  if(!is.null(Indels_tab_file)){
    if (file.exists(Indels_tab_file)){
      if(is.null(indels_obj)){
        indels_obj <- tabToIndelsClassification(indel.data = read.table(file = Indels_tab_file,sep = "\t",header = TRUE,
                                                                        check.names = FALSE,stringsAsFactors = FALSE),
                                                sampleID = sample_name,
                                                genome.v = genome.v)
      }else{
        message("[warning genomeChart] Indels_tab_file ignored because Indels already loaded from Indels_vcf_file")
        Indels_tab_file <- NULL
      }
    }else{
      message("[warning genomeChart] Indels_tab_file file not found: ",Indels_tab_file,". Ignoring and moving on.")
      Indels_tab_file <- NULL
    }
  }
  
  # Loading CNV file
  if(!is.null(CNV_tab_file)){
    if (file.exists(CNV_tab_file)){
      CNV_table <- read.table(file = CNV_tab_file,sep = "\t",header = TRUE,
                              check.names = FALSE,stringsAsFactors = FALSE)
    }else{
      message("[warning genomeChart] CNV_tab_file file not found: ",CNV_tab_file,". Ignoring and moving on.")
      CNV_tab_file <- NULL
    }
  }
  
  # Loading SV bedpe file
  if(!is.null(SV_bedpe_file)){
    if (file.exists(SV_bedpe_file)){
      sv_obj <- bedpeToRearrCatalogue(sv_bedpe = readTable(file = SV_bedpe_file))
      colnames(sv_obj$rearr_catalogue) <- "SV catalogue"
      if(!is.null(sv_obj$junctions_catalogue)) colnames(sv_obj$junctions_catalogue) <- "SV junction catalogue"
    }else{
      message("[warning genomeChart] SV_bedpe_file file not found: ",SV_bedpe_file,". Ignoring and moving on.")
      SV_bedpe_file <- NULL
    }
  }
  
  # check if there is anything at all to plot
  if(is.null(SNV_vcf_file) & is.null(SNV_tab_file) & is.null(Indels_vcf_file) & is.null(Indels_tab_file) & is.null(CNV_tab_file) & is.null(SV_bedpe_file)){
    message("[error genomeChart] no data found, nothing to plot. Quit.")
    return(NULL)
  }
  
  # get SBS catalogue if possible
  snvs_classified <- NULL
  if(!is.null(snvs_table)){
    sbs_obj <- tabToSNVcatalogue(subs = snvs_table,
                                 genome.v = genome.v)
    colnames(sbs_obj$muts)[c(1,2,4,5)] <- c("chr","position","REF","ALT")
    snvs_classified <- calcIntermutDist(sbs_obj$muts)
    colnames(sbs_obj$catalogue) <- "SNV catalogue"
  }
  
  # get DBS catalogue if possible
  if(!is.null(snvs_table)){
    dbs_obj <- tabToDNVcatalogue(muttable = snvs_table)
    colnames(dbs_obj$DNV_catalogue) <- "DNV catalogue"
  }
  
  # check SBS catalogues in SV clustering regions
  clustering_regions_sbs_catalogues <- NULL
  clusteringSBScatalogue_all <- NULL
  if(!is.null(sv_obj$clustering_regions) & !is.null(snvs_table)){
    # now find snvs in SV clusters and build catalogues
    resSBSinSVclusters <- getSBScataloguesInSVclusters(clustering_regions = sv_obj$clustering_regions,
                                                       snvs_table = snvs_table,
                                                       sample_name = sample_name,
                                                       genome.v = genome.v)
    clustering_regions_sbs_catalogues <- resSBSinSVclusters$clustering_regions_sbs_catalogues
    clusteringSBScatalogue_all <- resSBSinSVclusters$clusteringSBScatalogue_all
    colnames(clusteringSBScatalogue_all) <- "SNVs in SV clusters"
    snvs_table <- resSBSinSVclusters$snvs_table
  }
  
  # run the kataegis algorithm if requested
  kataegis_regions <- NULL
  kataegisSBScatalogue_all <- NULL
  if(!is.null(snvs_table) & runKataegis){
    kataegis <- findKataegis(snvs_table = snvs_table,
                             sample_name = sample_name)
    kataegis_regions <- kataegis$katregions
    snvs_table <- kataegis$snvs_table
    snvs_kataegis <- snvs_table[snvs_table$is.kataegis,,drop=F]
    if(nrow(snvs_kataegis)>0){
      resSBSkataegis <- tabToSNVcatalogue(subs = snvs_kataegis,
                                          genome.v = genome.v)
      kataegisSBScatalogue_all <- resSBSkataegis$catalogue
      colnames(kataegisSBScatalogue_all) <- "SNVs in kataegis"
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
  plotx <- 1.8
  ploty <- 1
  debug_grid_gap <- 0.1
  debug_grid_tick <- 0.01
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
    cairo_pdf(filename = outfilename,width = 7*plotx,height = 7*ploty)
  }else if(plottype=="png"){
    png(filename = outfilename,width = 2100*plotx,height = 2100*ploty,res = 300)
  }else{
    message("[error genomeChart] incorrect file type: ",plottype,". ",
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
  # par(fig = c(0,0.5,0,1),new=TRUE)
  placePanel(where = c(0,0),width = 1,height = 1)
  plotCircos(snvs_classified = snvs_classified,
             kataegis_regions = kataegis_regions,
             indels_table = indels_obj$indels_classified,
             indels_stats = indels_obj$count_proportion,
             CNV_table = CNV_table,
             sv_bedpe = sv_obj$annotated_bedpe,
             clustering_regions = sv_obj$clustering_regions,
             genome.v = genome.v)
  if(debug) drawDebugBox(1)
  
  # subs catalogue
  placePanel(where = c(0.96,0.73),width = 0.44,height = 0.21)
  if(!is.null(sbs_obj$catalogue)){
    plotSubsSignatures(sbs_obj$catalogue,
                       textscaling = 0.6)
  }else{
    plotMessage(msg = "SNV catalogue\nnot available",
                textscaling = msgtextscaling)
  }
  if(debug) drawDebugBox(2)
  
  # kataegis SBS catalogue
  placePanel(where = c(0.96,0.55),width = 0.44,height = 0.21)
  if(!is.null(kataegisSBScatalogue_all)){
    plotSubsSignatures(kataegisSBScatalogue_all,
                       textscaling = 0.6)
  }else{
    plotMessage(msg = "kataegis SBS catalogue\nnot available",
                textscaling = msgtextscaling)
  }
  if(debug) drawDebugBox(3)
  
  # SV clusters SBS catalogue
  placePanel(where = c(0.96,0.37),width = 0.44,height = 0.21)
  if(!is.null(clusteringSBScatalogue_all)){
    plotSubsSignatures(clusteringSBScatalogue_all,
                       textscaling = 0.6)
  }else{
    plotMessage(msg = "SV clusters SBS catalogue\nnot available",
                textscaling = msgtextscaling)
  }
  if(debug) drawDebugBox(4)
  
  # SV catalogue
  placePanel(where = c(1.35,0.68),width = 0.44,height = 0.26)
  if(!is.null(sv_obj$rearr_catalogue)){
    plotRearrSignatures(sv_obj$rearr_catalogue,
                        textscaling = 0.6,
                        mar = c(3.8, 3, 2, 1))
  }else{
    plotMessage(msg = "SV catalogue\nnot available",
                textscaling = msgtextscaling)
  }
  if(debug) drawDebugBox(5)
  
  # SV junctions catalogue
  placePanel(where = c(1.38,0.48),width = 0.38,height = 0.21)
  if(!is.null(sv_obj$junctions_catalogue)){
    plotJunctionsCatalogues(sv_obj$junctions_catalogue,
                            textscaling = 0.6,
                            mar = c(2, 3, 2, 1))
  }else{
    plotMessage(msg = "SV junction catalogue\nnot available",
                textscaling = msgtextscaling)
  }
  if(debug) drawDebugBox(6)
  
  # CN plot
  placePanel(where = c(1,0.01),width = 0.78,height = 0.21)
  if(!is.null(CNV_table)){
    plotCopyNumbers(sv_df = CNV_table,
                    sample_name = sample_name,
                    plottitle = "Copy Number Variations",
                    genome.v = genome.v,
                    mar = c(1.5,2.5,2,1),
                    textscaling = 0.6,
                    minorCNcolour = "#F08080FF",
                    totalCNcolour = "#9ACD32E6",
                    outofrangeCNcolour = "#604E97")
  }else{
    plotMessage(msg = "copy number data\nnot available",
                textscaling = msgtextscaling)
  }
  if(debug) drawDebugBox(7)
  
  # indels classification
  # placePanel(where = c(0.98,0.22),width = 0.4,height = 0.18)
  placePanel(where = c(1.39,0.3),width = 0.4,height = 0.18)
  if(!is.null(indels_obj$count_proportion)){
    plotIndelsClassSummary(indels_stats = indels_obj$count_proportion,
                           textscaling = 0.6,
                           mar = c(2,4,2,1))
  }else{
    plotMessage(msg = "Indels\nnot available",
                textscaling = msgtextscaling)
  }
  if(debug) drawDebugBox(8)
  
  # dnvs catalogue
  # placePanel(where = c(1.35,0.27),width = 0.44,height = 0.21)
  placePanel(where = c(0.95,0.19),width = 0.47,height = 0.21)
  if(!is.null(dbs_obj$DNV_catalogue)){
    plotDNVSignatures(convertToAlexandrovChannels(dbs_obj$DNV_catalogue),
                       textscaling = 0.6)
  }else{
    plotMessage(msg = "DNV catalogue\nnot available",
                textscaling = msgtextscaling)
  }
  if(debug) drawDebugBox(9)
  
  # placePanel(where = c(0.77,0.06),width = 0.21,height = 0.105)
  placePanel(where = c(1.55,0.21),width = 0.21,height = 0.105)
  plotCNVlegend(textscaling = 0.5)
  if(debug) drawDebugBox(10)
  
  # placePanel(where = c(0.08,0.06),width = 0.14,height = 0.1)
  placePanel(where = c(1.41,0.22),width = 0.14,height = 0.1)
  plotClustersLegend(textscaling = 0.5)
  if(debug) drawDebugBox(11)

  # close the file
  dev.off()
  
  # also I can return all the data/info obtained
  returnObj <- NULL
  returnObj$sample_name <- sample_name
  returnObj$genomev <- genome.v
  returnObj$snvs_table <- snvs_table
  returnObj$snvs_classified <- snvs_classified
  returnObj$sbs_obj <- sbs_obj
  returnObj$dbs_obj <- dbs_obj
  returnObj$kataegis_regions <- kataegis_regions
  returnObj$indels_obj <- indels_obj
  returnObj$CNV_table <- CNV_table
  returnObj$sv_obj <- sv_obj
  returnObj$clustering_regions_sbs_catalogues <- clustering_regions_sbs_catalogues
  returnObj$clusteringSBScatalogue_all <- clusteringSBScatalogue_all
  return(returnObj)
  
}

plotMessage <- function(msg,textscaling = 1){
  par(mar=c(0,0,0,0))
  plot(NULL,xlim = c(0,1),ylim = c(0,1),main="",type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",xaxs="i",yaxs="i")
  text(x = 0.5,
       y = 0.5,
       adj = 0.5,
       labels=msg,
       cex=textscaling)
}

getSBSCatalogueColours <- function(){
  snvcolours <- c(rgb(5, 195, 239, maxColorValue = 255),
                  rgb(0, 0, 0, maxColorValue = 255),
                  rgb(230, 47, 41, maxColorValue = 255),
                  rgb(208, 207, 207, maxColorValue = 255),
                  rgb(169, 212, 108, maxColorValue = 255),
                  rgb(238, 205, 204, maxColorValue = 255))
  names(snvcolours) <- c("C>A","C>G","C>T","T>A","T>C","T>G")
  return(snvcolours)
}

getIndelsClassColours <- function(){
  indels_colours <- c('firebrick4','firebrick1','firebrick3','darkgreen','grey')
  names(indels_colours) <- c("del.mhomology",
                             "del.repeatmediated",
                             "del.other",
                             "insertion",
                             "indel.complex")
  return(indels_colours)
}

getCNColoursAndRanges <- function(){
  CNcolours <- list()
  CNcolours$coloursTotalCN <- c("#D3D3D31A","#B3EE3A4C","#B3EE3A80","#B3EE3AB2",
                                "#B3EE3ABF","#9ACD32E6","#698B22E6")
  CNcolours$boundariesTotalCN <- c(0,3,4,8,16,32,64,10000)
  CNcolours$coloursMinorCN <- c("#F08080FF","#D3D3D31A")
  CNcolours$boundariesMinorCN <- c(0,1,10000)
  return(CNcolours)
}

getSVClassColours <- function(){
  sv_colours <- c("#1C86EEFF","#EE6A50FF","#006400FF","#595959FF")
  names(sv_colours) <- c("inversion",
                         "deletion",
                         "tandem-duplication",
                         "translocation")
  return(sv_colours)
}

getKellyColours <- function(){
  return(c('#F3C300', '#875692', '#F38400', '#A1CAF1', '#BE0032',
                    '#C2B280', '#848482', '#008856', '#E68FAC', '#0067A5', 
                    '#F99379', '#604E97', '#F6A600', '#B3446C', '#DCD300',
                    '#882D17', '#8DB600', '#654522', '#E25822',
                    '#2B3D26', '#222222', '#F2F3F4', '#CCCCCC'))
}

getClustersColours <- function(){
  kelly_colors <- getKellyColours()
  kataegis_region_colour <- kelly_colors[7]
  SVcluster_region_colour <- kelly_colors[7]
  svcols <- list()
  svcols$kataegis_region_colour <- kataegis_region_colour
  svcols$SVcluster_region_colour <- SVcluster_region_colour
  return(svcols)
}

plotIndelsClassSummary <- function(indels_stats,
                                   textscaling = 1,
                                   mar = c(3,6,2,2)){
  indels_colours <- getIndelsClassColours()
  indels_colours <- indels_colours[c("del.mhomology","del.repeatmediated","del.other","insertion","indel.complex")]
  values_to_plot <- indels_stats[,c("del.mh","del.rep","del.none","ins","complex")]
  plot_labels <- c("deletion at MH","deletion at repeat","deletion other","insertion","complex")
  par(mar=mar)
  barplot(t(as.matrix(values_to_plot)),beside = T,
          horiz = T,space = 0.2,
          col = indels_colours,
          names.arg = plot_labels,
          mar=mar,
          las=2,
          main = paste0(" insertions and deletions\n(",indels_stats$all.indels," Indels)"),
          cex.axis = textscaling*0.7,
          cex.names = textscaling*0.7,
          border = NA,
          cex.main=textscaling)
}

plotCNVlegend <- function(textscaling = 1){
  # get colours and ranges
  CNVcolours <- getCNColoursAndRanges()
  
  # plot parameters
  leftmarginblocks <- 1
  xgapmarginblocks <- 3
  rightmarginblocks <- 1
  topmarginblocks <- 3
  bottommarginblocks <- 3
  nxblocks <- length(CNVcolours$coloursTotalCN) + length(CNVcolours$coloursMinorCN) + leftmarginblocks + xgapmarginblocks + rightmarginblocks
  nyblocks <- topmarginblocks + bottommarginblocks + 1
  
  # plot it
  par(mar=c(0,0,0,0))
  plot(NULL,xlim = c(0,nxblocks),ylim = c(0,nyblocks),main="",type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",xaxs="i",yaxs="i")
  
  # total CN legend (gain)
  text(x = leftmarginblocks + length(CNVcolours$coloursTotalCN)/2,
       y = bottommarginblocks + 2.5,
       labels = "\nTotal CN",
       cex = textscaling)
  
  for(i in 1:length(CNVcolours$coloursTotalCN)){
    rect(xleft = leftmarginblocks + (i-1),
         ybottom = bottommarginblocks,
         xright = leftmarginblocks + i,
         ytop = bottommarginblocks+1,
         col = CNVcolours$coloursTotalCN[i],
         border = NA)
    if(i < length(CNVcolours$coloursTotalCN)){
      cnlow <- CNVcolours$boundariesTotalCN[i]
      cnhigh <- CNVcolours$boundariesTotalCN[i+1]-1
      if(cnlow < cnhigh){
        textlabel <- paste0(cnlow," - ",cnhigh)
      }else if(cnlow==cnhigh){
        textlabel <- cnlow
      }
    }else{
      textlabel <- paste0("\u2265 ",CNVcolours$boundariesTotalCN[i])
    }
    
    text(x = leftmarginblocks + (i-1) + 0.5,
         y = bottommarginblocks - 0.5,srt = 90,adj = 1,
         labels = textlabel,
         cex = textscaling)
  }
  
  # minor CN legend (LOH)
  leftstart <- leftmarginblocks + length(CNVcolours$coloursTotalCN) + xgapmarginblocks
  
  text(x = leftstart + length(CNVcolours$coloursMinorCN)/2,
       y = bottommarginblocks + 2.5,
       labels = "\nMinor CN",
       cex = textscaling)
  
  for(i in 1:length(CNVcolours$coloursMinorCN)){
    rect(xleft = leftstart + (i-1),
         ybottom = bottommarginblocks,
         xright = leftstart + i,
         ytop = bottommarginblocks+1,
         col = CNVcolours$coloursMinorCN[i],
         border = NA)
    if(i < length(CNVcolours$coloursMinorCN)){
      cnlow <- CNVcolours$boundariesMinorCN[i]
      cnhigh <- CNVcolours$boundariesMinorCN[i+1]-1
      if(cnlow < cnhigh){
        textlabel <- paste0(cnlow," - ",cnhigh)
      }else if(cnlow==cnhigh){
        textlabel <- cnlow
      }
    }else{
      textlabel <- paste0("\u2265 ",CNVcolours$boundariesMinorCN[i])
    }
    
    text(x = leftstart + (i-1) + 0.5,
         y = bottommarginblocks - 0.5,srt = 90,adj = 1,
         labels = textlabel,
         cex = textscaling)
  }
}

plotClustersLegend <- function(textscaling=1){
  svcols <- getClustersColours()
  kataegis_region_colour <- svcols$kataegis_region_colour
  SVcluster_region_colour <- svcols$SVcluster_region_colour
  # plot it
  par(mar=c(0,0,0,0))
  plot(NULL,xlim = c(0,7),ylim = c(0,5),main="",type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",xaxs="i",yaxs="i")
  # kataegis
  lines(x=c(0.5,2.5),
        y=c(3,3),
        col=svcols$kataegis_region_colour,
        lwd=1)
  lines(x=c(0.5,1.5),
        y=c(3,3),
        col=svcols$kataegis_region_colour,
        lwd=2)
  text(x=3,y=3,labels="kataegis",cex=textscaling,adj=0)
  # svclusters
  lines(x=c(2,1,1,2),
        y=c(1.1,1.1,1.9,1.9),
        col=svcols$SVcluster_region_colour,
        lwd=2)
  text(x=3,y=1.5,labels="SV clusters",cex=textscaling,adj=0)
}


plotCircos <- function(snvs_classified,
                       kataegis_regions,
                       indels_table,
                       indels_stats,
                       CNV_table,
                       sv_bedpe,
                       clustering_regions,
                       genome.v){
  
  # set some colours
  kelly_colors <- getKellyColours()
  # SNVs colours
  snvcolours <- getSBSCatalogueColours()
  
  # Indel colours
  indels_colours <- getIndelsClassColours()
  
  # CNV colours
  CNcolours <- getCNColoursAndRanges()
  coloursTotalCN <- CNcolours$coloursTotalCN
  boundariesTotalCN <- CNcolours$boundariesTotalCN
  coloursMinorCN <- CNcolours$coloursMinorCN
  boundariesMinorCN <- CNcolours$boundariesMinorCN
  
  getCNVcolour <- function(CNVval,coloursScale,boundariesScale){
    found <- FALSE
    i <- 1
    while (!found) {
      if(CNVval >= boundariesScale[i] & CNVval < boundariesScale[i+1]){
        found <- TRUE
      }else{
        i <- i + 1
      }
    }
    return(coloursScale[i])
  }
  
  # svclass colours
  sv_colours <- getSVClassColours()
  
  # clusters colours
  svcols <- getClustersColours()
  kataegis_region_colour <- svcols$kataegis_region_colour
  SVcluster_region_colour <- svcols$SVcluster_region_colour
                             
  if(!is.null(snvs_classified)){
    if(nrow(snvs_classified)>0){
      if(!startsWith(as.character(snvs_classified$chr[1]),"chr")){
        snvs_classified$chr <- paste0("chr",snvs_classified$chr)
      }
    }
  }
  
  if(!is.null(kataegis_regions)){
    if(nrow(kataegis_regions)>0){
      if(!startsWith(as.character(kataegis_regions$chr[1]),"chr")){
        kataegis_regions$chr <- paste0("chr",kataegis_regions$chr)
      }
    }
  }
  
  if(!is.null(indels_table)){
    if(nrow(indels_table)>0){
      if(!startsWith(as.character(indels_table$chr[1]),"chr")){
        indels_table$chr <- paste0("chr",indels_table$chr)
      }
    }
  }
  
  if(!is.null(CNV_table)){
    if(nrow(CNV_table)>0){
      if(!startsWith(as.character(CNV_table$Chromosome[1]),"chr")){
        CNV_table$Chromosome <- paste0("chr",CNV_table$Chromosome)
      }
    }
  }
  
  if(!is.null(sv_bedpe)){
    if(nrow(sv_bedpe)>0){
      if(!startsWith(as.character(sv_bedpe$chrom1[1]),"chr")){
        sv_bedpe$chrom1 <- paste0("chr",sv_bedpe$chrom1)
        sv_bedpe$chrom2 <- paste0("chr",sv_bedpe$chrom2)
      }
    }
  }
  
  if(!is.null(clustering_regions)) {
    if(!startsWith(as.character(clustering_regions$chr[1]),"chr")){
      clustering_regions$chr <- paste0("chr",clustering_regions$chr)
    }
  }
  
  par(mai=c(0,0,0,0))
  circlize::circos.par("start.degree" = 90)
  circlize::circos.initializeWithIdeogram(species = genome.v,plotType = NULL)
  
  # chromosome names
  circlize::circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chromNumber <- gsub("chr", "", circlize::CELL_META$sector.index)
    circlize::circos.text(circlize::CELL_META$xcenter, circlize::CELL_META$ylim[1]+(circlize::CELL_META$ylim[2]-circlize::CELL_META$ylim[1])*0.3, 
                          chromNumber, cex = 0.6, niceFacing = TRUE)
  }, track.height = 0.05, cell.padding = c(0, 0, 0, 0), bg.border = NA)
  
  circlize::circos.genomicIdeogram(species = genome.v)
  
  # draw rectangles where the kataegis clusters are
  circlize::circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chromNumber <- circlize::CELL_META$sector.index
    
    if(!is.null(kataegis_regions)){
      tmpRows <- kataegis_regions[kataegis_regions$chr==chromNumber,,drop=F]
      if(nrow(tmpRows)>0){
        for (i in 1:nrow(tmpRows)){
          circlize::circos.rect(xleft = tmpRows$start.bp[i],
                                ybottom = circlize::CELL_META$ylim[1],
                                xright = tmpRows$end.bp[i],
                                ytop = circlize::CELL_META$ylim[2],
                                border = kataegis_region_colour,
                                col=kataegis_region_colour,lwd = 2)
        }
      }
    }
  }, track.height = 0.02, cell.padding = c(0, 0, 0, 0), bg.border = NA,
     bg.lwd = 0.5, track.margin = c(0.005,0.01))
  
  # draw snvs
  maxLogDistance <- 8
  circlize::circos.track(ylim = c(0, maxLogDistance), panel.fun = function(x, y) {
    chromNumber <- circlize::CELL_META$sector.index
    
    if(!is.null(kataegis_regions)){
      tmpRows <- kataegis_regions[kataegis_regions$chr==chromNumber,,drop=F]
      if(nrow(tmpRows)>0){
        for (i in 1:nrow(tmpRows)){
          meanpos <- tmpRows$start.bp[i] + (tmpRows$end.bp[i] - tmpRows$start.bp[i])/2
          suppressMessages(circlize::circos.lines(x = c(meanpos,meanpos),
                                                  y = c(0,maxLogDistance+1),
                                                  col = kataegis_region_colour))
        }
      }
    }
    
    if(!is.null(snvs_classified)){
      tmpRows <- snvs_classified[snvs_classified$chr==chromNumber,,drop=F]
      if(nrow(tmpRows)>0){
        pos <- tmpRows$position
        logdist <- log10(tmpRows$distPrev)
        snvtype <- paste(tmpRows$pyrwt,tmpRows$pyrmut,sep = ">")
        circlize::circos.points(pos, logdist,pch = 16,cex = 0.4,col = snvcolours[snvtype])
      }
    }
  }, track.height = 0.2, cell.padding = c(0, 0, 0, 0), bg.border = "lightgrey",
     bg.lwd = 0.5, track.margin = c(0.005,0))
  
  # draw insertions
  circlize::circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chromNumber <- circlize::CELL_META$sector.index
    if(!is.null(indels_table)){
      tmpRows <- indels_table[indels_table$chr==chromNumber & indels_table$indel.type=="I",,drop=F]
      if(nrow(tmpRows)>0){
        for(i in 1:nrow(tmpRows)){
          circlize::circos.lines(x = c(tmpRows$pos[i],tmpRows$pos[i]),
                                 y = c(0,1),
                                 col = indels_colours[tmpRows$indel.class[i]])
        }
      }
    }
  }, track.height = 0.04, cell.padding = c(0, 0, 0, 0), bg.border = "lightgrey",
  bg.lwd = 0.5, track.margin = c(0.005,0.005))
  
  # draw deletions
  circlize::circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chromNumber <- circlize::CELL_META$sector.index
    if(!is.null(indels_table)){
      tmpRows <- indels_table[indels_table$chr==chromNumber & indels_table$indel.type=="D",,drop=F]
      if(nrow(tmpRows)>0){
        for(i in 1:nrow(tmpRows)){
          circlize::circos.lines(x = c(tmpRows$pos[i],tmpRows$pos[i]),
                                 y = c(0,1),
                                 col = indels_colours[tmpRows$indel.class[i]])
        }
      }
    }
  }, track.height = 0.04, cell.padding = c(0, 0, 0, 0), bg.border = "lightgrey",
  bg.lwd = 0.5, track.margin = c(0.005,0.005))
  
  # draw Total copy number data (gains)
  circlize::circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chromNumber <- circlize::CELL_META$sector.index
    
    if(!is.null(CNV_table)){
      tmpRows <- CNV_table[CNV_table$Chromosome==chromNumber,,drop=F]
      if(nrow(tmpRows)>0){
        for (i in 1:nrow(tmpRows)){
          circlize::circos.rect(xleft = tmpRows$chromStart[i],
                                ybottom = circlize::CELL_META$ylim[1],
                                xright = tmpRows$chromEnd[i],
                                ytop = circlize::CELL_META$ylim[2],
                                border = NA,
                                col=getCNVcolour(CNVval = tmpRows$total.copy.number.inTumour[i],
                                                 coloursScale = coloursTotalCN,
                                                 boundariesScale = boundariesTotalCN),
                                lwd=1)
        }
      }
    }
  }, track.height = 0.04, cell.padding = c(0, 0, 0, 0), bg.border = NA,
  bg.lwd = 0.5, track.margin = c(0.005,0.005))
  
  # draw Minor copy number data (LOH)
  circlize::circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chromNumber <- circlize::CELL_META$sector.index
    
    if(!is.null(CNV_table)){
      tmpRows <- CNV_table[CNV_table$Chromosome==chromNumber,,drop=F]
      if(nrow(tmpRows)>0){
        for (i in 1:nrow(tmpRows)){
          circlize::circos.rect(xleft = tmpRows$chromStart[i],
                                ybottom = circlize::CELL_META$ylim[1],
                                xright = tmpRows$chromEnd[i],
                                ytop = circlize::CELL_META$ylim[2],
                                border = NA,
                                col=getCNVcolour(CNVval = tmpRows$minor.copy.number.inTumour[i],
                                                 coloursScale = coloursMinorCN,
                                                 boundariesScale = boundariesMinorCN),
                                lwd=1)
        }
      }
    }
  }, track.height = 0.04, cell.padding = c(0, 0, 0, 0), bg.border = NA,
  bg.lwd = 0.5, track.margin = c(0.005,0.005))
  
  # draw rectangles where the SV clusters are
  circlize::circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chromNumber <- circlize::CELL_META$sector.index
    
    if(!is.null(clustering_regions)){
      tmpRows <- clustering_regions[clustering_regions$chr==chromNumber,,drop=F]
      if(nrow(tmpRows)>0){
        for (i in 1:nrow(tmpRows)){
          circlize::circos.lines(x = c(tmpRows$start.bp[i],tmpRows$start.bp[i],
                                       tmpRows$end.bp[i],tmpRows$end.bp[i]),
                                 y = c(0,1,1,0),
                                 col = SVcluster_region_colour)
        }
      }
    }
  }, track.height = 0.02, cell.padding = c(0, 0, 0, 0), bg.border = NA,
     bg.lwd = 0.5, track.margin = c(0,0.005))
  
  # draw SV links
  if(!is.null(sv_bedpe)){
    for(j in 1:nrow(sv_bedpe)){
      circlize::circos.link(sector.index1 = as.character(sv_bedpe$chrom1[j]),
                            point1 = sv_bedpe$start1[j],
                            sector.index2 = as.character(sv_bedpe$chrom2[j]),
                            point2 = sv_bedpe$start2[j],
                            col = sv_colours[sv_bedpe$svclass[j]],
                            rou = circlize::get_most_inside_radius() + 0.005)
      
    }
  }
 
  circlize::circos.clear()
  
}

