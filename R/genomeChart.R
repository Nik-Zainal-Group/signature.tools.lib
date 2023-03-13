
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
                        genome.v = "hg19"){
  
  # set up some variables
  snvs_table <- NULL
  sbs_obj <- NULL
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
    snvs_table <- resSBSinSVclusters$snvs_table
  }
  
  # run the kataegis algorithm if requested
  kataegis_regions <- NULL
  if(!is.null(snvs_table) & runKataegis){
    kataegis <- findKataegis(snvs_table = snvs_table,
                             sample_name = sample_name)
    kataegis_regions <- kataegis$katregions
    snvs_table <- kataegis$snvs_table
  }
  
  # now plot
  plotCircos(snvs_classified = snvs_classified,
             kataegis_regions = kataegis_regions,
             indels_table = indels_obj$indels_classified,
             indels_stats = indels_obj$count_proportion,
             CNV_table = CNV_table,
             sv_bedpe = sv_obj$annotated_bedpe,
             clustering_regions = sv_obj$clustering_regions,
             genome.v = genome.v)
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
  kelly_colors <- c('#F3C300', '#875692', '#F38400', '#A1CAF1', '#BE0032',
                    '#C2B280', '#848482', '#008856', '#E68FAC', '#0067A5', 
                    '#F99379', '#604E97', '#F6A600', '#B3446C', '#DCD300',
                    '#882D17', '#8DB600', '#654522', '#E25822',
                    '#2B3D26', '#222222', '#F2F3F4', '#CCCCCC')
  # SNVs colours
  snvcolours <- c(rgb(5, 195, 239, maxColorValue = 255),
                  rgb(0, 0, 0, maxColorValue = 255),
                  rgb(230, 47, 41, maxColorValue = 255),
                  rgb(208, 207, 207, maxColorValue = 255),
                  rgb(169, 212, 108, maxColorValue = 255),
                  rgb(238, 205, 204, maxColorValue = 255))
  names(snvcolours) <- c("C>A","C>G","C>T","T>A","T>C","T>G")
  
  # Indel colours
  indels_colours <- c('firebrick4','firebrick1','firebrick3','darkgreen','grey')
  names(indels_colours) <- c("del.mhomology",
                             "del.repeatmediated",
                             "del.other",
                             "insertion",
                             "indel.complex")
  
  # CNV colours
  coloursTotalCN <- c( "#D3D3D31A","#B3EE3A4C","#B3EE3A80","#B3EE3AB2",
                       "#B3EE3ABF","#9ACD32E6","#698B22E6")
  boundariesTotalCN <- c(0,3,4,8,16,32,64,10000)
  coloursMinorCN <- c("#F08080FF","#D3D3D31A")
  boundariesMinorCN <- c(0,1,10000)
  
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
  sv_colours <- c("#1C86EEFF","#EE6A50FF","#006400FF","#595959FF")
  names(sv_colours) <- c("inversion",
                         "deletion",
                         "tandem-duplication",
                         "translocation")
  
  # clusters colours
  kataegis_region_colour <- kelly_colors[21]
  SVcluster_region_colour <- kelly_colors[21]
                             
  if(!is.null(snvs_classified)){
    if(nrow(snvs_classified)>0){
      if(!startsWith(snvs_classified$chr[1],"chr")){
        snvs_classified$chr <- paste0("chr",snvs_classified$chr)
      }
    }
  }
  
  if(!is.null(kataegis_regions)){
    if(nrow(kataegis_regions)>0){
      if(!startsWith(kataegis_regions$chr[1],"chr")){
        kataegis_regions$chr <- paste0("chr",kataegis_regions$chr)
      }
    }
  }
  
  if(!is.null(indels_table)){
    if(nrow(indels_table)>0){
      if(!startsWith(indels_table$chr[1],"chr")){
        indels_table$chr <- paste0("chr",indels_table$chr)
      }
    }
  }
  
  if(!is.null(CNV_table)){
    if(nrow(CNV_table)>0){
      if(!startsWith(CNV_table$Chromosome[1],"chr")){
        CNV_table$Chromosome <- paste0("chr",CNV_table$Chromosome)
      }
    }
  }
  
  if(!is.null(sv_bedpe)){
    if(nrow(sv_bedpe)>0){
      if(!startsWith(sv_bedpe$chrom1[1],"chr")){
        sv_bedpe$chrom1 <- paste0("chr",sv_bedpe$chrom1)
        sv_bedpe$chrom2 <- paste0("chr",sv_bedpe$chrom2)
      }
    }
  }
  
  if(!is.null(clustering_regions)) {
    if(!startsWith(clustering_regions$chr[1],"chr")){
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
  }, track.height = 0.01, cell.padding = c(0, 0, 0, 0), bg.border = NA,
     bg.lwd = 0.5, track.margin = c(0.005,0.01))
  
  # draw snvs
  circlize::circos.track(ylim = c(0, 8), panel.fun = function(x, y) {
    chromNumber <- circlize::CELL_META$sector.index
    if(!is.null(snvs_classified)){
      tmpRows <- snvs_classified[snvs_classified$chr==chromNumber,,drop=F]
      if(nrow(tmpRows)>0){
        pos <- tmpRows$position
        logdist <- log10(tmpRows$distPrev)
        snvtype <- paste(tmpRows$pyrwt,tmpRows$pyrmut,sep = ">")
        circlize::circos.points(pos, logdist,pch = 16,cex = 0.3,col = snvcolours[snvtype])
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
  }, track.height = 0.05, cell.padding = c(0, 0, 0, 0), bg.border = "lightgrey",
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
  }, track.height = 0.05, cell.padding = c(0, 0, 0, 0), bg.border = "lightgrey",
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
  }, track.height = 0.05, cell.padding = c(0, 0, 0, 0), bg.border = NA,
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
  }, track.height = 0.05, cell.padding = c(0, 0, 0, 0), bg.border = NA,
  bg.lwd = 0.5, track.margin = c(0.005,0.005))
  
  # draw rectangles where the SV clusters are
  circlize::circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chromNumber <- circlize::CELL_META$sector.index
    
    if(!is.null(clustering_regions)){
      tmpRows <- clustering_regions[clustering_regions$chr==chromNumber,,drop=F]
      if(nrow(tmpRows)>0){
        for (i in 1:nrow(tmpRows)){
          circlize::circos.rect(xleft = tmpRows$start.bp[i],
                                ybottom = circlize::CELL_META$ylim[1],
                                xright = tmpRows$end.bp[i],
                                ytop = circlize::CELL_META$ylim[2],
                                border = SVcluster_region_colour,
                                col=SVcluster_region_colour,
                                lwd=2)
        }
      }
    }
  }, track.height = 0.01, cell.padding = c(0, 0, 0, 0), bg.border = NA,
     bg.lwd = 0.5, track.margin = c(0,0.005))
  
  # draw SV links
  if(!is.null(sv_bedpe)){
    for(j in 1:nrow(sv_bedpe)){
      circlize::circos.link(sector.index1 = sv_bedpe$chrom1[j],
                            point1 = sv_bedpe$start1[j],
                            sector.index2 = sv_bedpe$chrom2[j],
                            point2 = sv_bedpe$start2[j],
                            col = sv_colours[sv_bedpe$svclass[j]],
                            rou = circlize::get_most_inside_radius() + 0.005)
      
    }
  }
 
  circlize::circos.clear()
}

