#' genomeChartSV
#' 
#' Plotting of somatic variants using the circlize R package. This plot focuses
#' on structural variants (SVs or rearrangements). Clusters of SVs, if any, are
#' shown on the plot using a different colour for each cluster. If SNVs are
#' available, SNVs within the SV clusters are identified and used to build an
#' SNV mutational catalogue of these localised variants. Moreover, if the bedpe
#' file of the SVs contains data about non-templated insertions and micro-homology
#' deletions at breakpoint junctions then the SV junction catalogue is also
#' built and shown.
#' 
#' @param outfilename file where the figure should be plotted. This also determines the file type of the output, use either pdf (recommended) or png
#' @param SV_bedpe_file name of the tab separated bedpe file containing the SVs. The file should contain a rearrangement for each row (two breakpoint positions should be on one row as determined by a pair of mates of paired-end sequencing) and should already be filtered according to the user preference, as all rearrangements in the file will be used and no filter will be applied. The files should contain a header in the first line with the following columns: "chrom1", "start1", "end1", "chrom2", "start2", "end2" and "sample" (sample name). In addition, either two columns indicating the strands of the mates, "strand1" (+ or -) and "strand2" (+ or -), or one column indicating the structural variant class, "svclass": translocation, inversion, deletion, tandem-duplication. The column "svclass" should correspond to (Sanger BRASS convention): inversion (strands +/- or -/+ and mates on the same chromosome), deletion (strands +/+ and mates on the same chromosome), tandem-duplication (strands -/- and mates on the same chromosome), translocation (mates are on different chromosomes). In addition, columns 'non-template'	and 'micro-homology' can be specified, including the non-templated insertion or micro-homology deletion sequence, which will be used to build an SV junctions catalogue.
#' @param SNV_vcf_file name of the vcf file containing the SNVs
#' @param SNV_tab_file name of the tab separated file containing the SNVs. Column names should be: chr, position, REF and ALT. If SNV_vcf_file is also specified, these variants will be ignored and the variants in the vcf file will be used instead.
#' @param plot_title title of the plot. If NULL, then the sample_name will be used as title. Use "", the empty string, to have no title.
#' @param genome.v genome version to use, either hg19 or hg38
#' @return all computed results, like catalogues and clustering regions, will be returned
#' @export
genomeChartSV <- function(outfilename,
                          SV_bedpe_file,
                          SNV_vcf_file = NULL,
                          SNV_tab_file = NULL,
                          PEAK.FACTOR = 10,
                          kmin = 10,
                          plot_title = NULL,
                          genome.v = "hg19"){
  # some checks first
  if (!file.exists(SV_bedpe_file)){
    message("[error genomeChartSV] SV_bedpe_file file not found: ",SV_bedpe_file)
    return(NULL)
  }
  # load the bedpe
  sv_bedpe <- readTable(file = SV_bedpe_file)
  
  # annotate structural variants and retrieve the rearrangement catalogue
  # this will also check that the table is in the correct format
  sv_obj <- bedpeToRearrCatalogue(sv_bedpe = sv_bedpe,
                                  kmin = kmin,
                                  PEAK.FACTOR = PEAK.FACTOR)
  
  sample_name <- unique(sv_obj$annotated_bedpe$sample)
  
  # check if we have SNVs
  snvs_table <- NULL
  if(!is.null(SNV_vcf_file)){
    if (file.exists(SNV_vcf_file)){
      snvs_table <- fromVcfToTable(vcfFilename = SNV_vcf_file,
                                   genome.v = genome.v)
    }else{
      message("[warning genomeChartSV] SNV_vcf_file file not found: ",SNV_vcf_file,". Ignoring and moving on.")
    }
  }else if(!is.null(SNV_tab_file)){
    if (file.exists(SNV_tab_file)){
      snvs_table <- read.table(file = SNV_tab_file,sep = "\t",header = TRUE,
                               check.names = FALSE,stringsAsFactors = FALSE)
    }else{
      message("[warning genomeChartSV] SNV_tab_file file not found: ",SNV_tab_file,". Ignoring and moving on.")
    }
  }

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
  
  # time to plot, outfilename needs to be specified
  plottype <- substr(outfilename, nchar(outfilename) - 2, nchar(outfilename))
  dir.create(dirname(outfilename),showWarnings = F,recursive = T)
  
  # open the file
  if(plottype=="pdf"){
    cairo_pdf(filename = outfilename,width = 10,height = 7)
  }else if(plottype=="png"){
    png(filename = outfilename,width = 3000,height = 2100,res = 300)
  }else{
    message("[error genomeChartSV] incorrect file type: ",plottype,". ",
            "Please end your file name with .pdf or .png")
    return(NULL)
  }

  # plot data
  if(is.null(plot_title)) plot_title <- sample_name
  par(fig = c(0,1,0,1),mai=c(0,0,0.5,0))
  plot(NULL,xlim = c(0,1),ylim = c(0,1),main = plot_title,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
  par(fig = c(0,0.7,0,1),new=TRUE)
  plotCircosSV(sv_bedpe = sv_obj$annotated_bedpe,
               clustering_regions = sv_obj$clustering_regions,
               genome.v = genome.v)
  par(fig = c(0.65,0.98,0.5,0.9),new=TRUE)
  plotRearrSignatures(sv_obj$rearr_catalogue,textscaling = 0.6)
  par(fig = c(0.69,0.97,0.23,0.63),new=TRUE)
  plotJunctionsCatalogues(sv_obj$junctions_catalogue,textscaling = 0.6)
  if(!is.null(clusteringSBScatalogue_all)){
    par(fig = c(0.65,0.98,0.17,0.37),new=TRUE)
    plotSubsSignatures(clusteringSBScatalogue_all,textscaling = 0.6)
  }
  # close the file
  dev.off()
  
  # also I can return all the data/info obtained
  returnObj <- NULL
  returnObj$sv_obj <- sv_obj
  returnObj$clustering_regions_sbs_catalogues <- clustering_regions_sbs_catalogues
  returnObj$clusteringSBScatalogue_all <- clusteringSBScatalogue_all
  returnObj$snvs_table <- snvs_table
  return(returnObj)
}

fromVcfToTable <- function(vcfFilename,
                           genome.v){
  vcf_data <- VariantAnnotation::readVcf(vcfFilename, genome.v)
  vcf_data <- VariantAnnotation::expand(vcf_data)
  rd <- SummarizedExperiment::rowRanges(vcf_data)
  ref <- as.character(rd$REF)
  alt <- as.character(rd$ALT)
  chr <- as.character(GenomeInfoDb::seqnames(vcf_data))
  rgs <- IRanges::ranges(vcf_data)
  position <- BiocGenerics::start(rgs)
  muts_table <- data.frame(chr = chr,
                           position = position,
                           REF = ref,
                           ALT = alt,
                           stringsAsFactors = F)
  return(muts_table)
}


getSVinRange <- function(sv_bedpe,chr,pstart,pend){
  bp1in <- sv_bedpe$chrom1==chr & sv_bedpe$start1>=pstart & sv_bedpe$start1<=pend
  bp2in <- sv_bedpe$chrom2==chr & sv_bedpe$start2>=pstart & sv_bedpe$start2<=pend
  return(bp1in | bp2in)
}

getSBScataloguesInSVclusters <- function(clustering_regions,
                                         snvs_table,
                                         sample_name,
                                         genome.v){
  
  clustering_regions_sbs_catalogues <- list()
  snvs_table$SVcluster <- NA
  
  for(i in 1:nrow(clustering_regions)){
    region_chr <- clustering_regions$chr[i]
    region_start <- clustering_regions$start.bp[i]
    region_end <- clustering_regions$end.bp[i]
    select_snvs <- snvs_table$chr==region_chr & snvs_table$position >= region_start & snvs_table$position <= region_end
    region_snvs <- snvs_table[select_snvs,]
    snvs_table[select_snvs,"SVcluster"] <- i
    res_sbs <- tabToSNVcatalogue(region_snvs,
                                 genome.v = genome.v)
    colnames(res_sbs$catalogue) <- paste0(sample_name," - SV cluster ",i)
    clustering_regions_sbs_catalogues[[i]] <- res_sbs$catalogue
  }
  clustering_regions_sbs_catalogues <- do.call(cbind,clustering_regions_sbs_catalogues)
  clusteringSBScatalogue_all <- data.frame(apply(clustering_regions_sbs_catalogues, 1, sum),stringsAsFactors = F)
  colnames(clusteringSBScatalogue_all) <- paste0(sample_name," - SNVs in SV clusters")
  
  res <- list()
  res$clusteringSBScatalogue_all <- clusteringSBScatalogue_all
  res$clustering_regions_sbs_catalogues <- clustering_regions_sbs_catalogues
  res$snvs_table <- snvs_table
  return(res)
}

plotCircosSV <- function(sv_bedpe,
                         clustering_regions,
                         genome.v){
  kelly_colors <- c('#F3C300', '#875692', '#F38400', '#A1CAF1', '#BE0032',
                    '#C2B280', '#848482', '#008856', '#E68FAC', '#0067A5', 
                    '#F99379', '#604E97', '#F6A600', '#B3446C', '#DCD300',
                    '#882D17', '#8DB600', '#654522', '#E25822',
                    '#2B3D26', '#222222', '#F2F3F4', '#CCCCCC')
  
  if(nrow(sv_bedpe)>0){
    if(!startsWith(sv_bedpe$chrom1[1],"chr")){
      sv_bedpe$chrom1 <- paste0("chr",sv_bedpe$chrom1)
      sv_bedpe$chrom2 <- paste0("chr",sv_bedpe$chrom2)
    }
  }
  
  if(!is.null(clustering_regions)) {
    maxcolours <- length(kelly_colors)
    if(nrow(clustering_regions)>maxcolours){
      clustering_regions$colour <- rep("#CCCCCC",nrow(clustering_regions))
      clustering_regions$colour[1:maxcolours] <- kelly_colors[1:maxcolours]
    }else{
      clustering_regions$colour <- kelly_colors[1:nrow(clustering_regions)]
    }
    if(!startsWith(clustering_regions$chr[1],"chr")){
      clustering_regions$chr <- paste0("chr",clustering_regions$chr)
    }
  }
  
  circlize::circos.par("start.degree" = 90)
  circlize::circos.initializeWithIdeogram(species = genome.v,plotType = NULL)
  
  # chromosome names
  circlize::circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chromNumber <- gsub("chr", "", circlize::CELL_META$sector.index)
    circlize::circos.text(circlize::CELL_META$xcenter, circlize::CELL_META$ylim[1]+(circlize::CELL_META$ylim[2]-circlize::CELL_META$ylim[1])*0.3, 
                          chromNumber, cex = 0.6, niceFacing = TRUE)
  }, track.height = 0.05, cell.padding = c(0, 0, 0, 0), bg.border = NA)
  
  circlize::circos.genomicIdeogram(species = genome.v)
  
  # draw rectangles where the clusters are
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
                                border = NA,col=tmpRows$colour[i])
        }
      }
    }
  }, track.height = 0.1, cell.padding = c(0, 0, 0, 0), bg.border = "lightgrey", bg.lwd = 0.5)
  
  # draw unclustered links first
  tmpBEDPE <- sv_bedpe[sv_bedpe$is.clustered==F,,drop=F]
  for(j in 1:nrow(tmpBEDPE)){
    circlize::circos.link(sector.index1 = tmpBEDPE$chrom1[j],
                          point1 = tmpBEDPE$start1[j],
                          sector.index2 = tmpBEDPE$chrom2[j],
                          point2 = tmpBEDPE$start2[j],col = "lightgrey")

  }
  
  # draw links from clusters
  if(!is.null(clustering_regions)){
    for(i in 1:nrow(clustering_regions)){
      #i <- 1
      whichSVs <- getSVinRange(sv_bedpe,clustering_regions$chr[i],clustering_regions$start.bp[i],clustering_regions$end.bp[i])
      if(sum(whichSVs)>0){
        tmpBEDPE <- sv_bedpe[whichSVs,,drop=F]
        for(j in 1:nrow(tmpBEDPE)){
          circlize::circos.link(sector.index1 = tmpBEDPE$chrom1[j],
                                point1 = tmpBEDPE$start1[j],
                                sector.index2 = tmpBEDPE$chrom2[j],
                                point2 = tmpBEDPE$start2[j],col = clustering_regions$colour[i],lwd = 2)
        }
        
      }
      
    }
  }
  circlize::circos.clear()
}


