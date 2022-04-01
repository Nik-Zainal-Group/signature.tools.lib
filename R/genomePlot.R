#  plotting colours and parameters
set.plot.params <- function(colour.scheme = "ascat"){

  # circos parameters
  params.my <- list()

  #use these two params to adjust circle size
  params.my$plot.radius <- 2.15
  params.my$genomeplot.margin <- 0.25

  params.my$track.background <- 'white'
  params.my$highlight.width <- 0.2
  params.my$point.size <- 0.3
  params.my$point.type <- 16
  params.my$radius.len <- 3
  params.my$chr.ideog.pos <- 3.2
  params.my$highlight.pos <- 2.09 #3.35
  params.my$chr.name.pos <- 2.14 #3.45
  params.my$track.in.start <- 3.05
  params.my$track.out.start <- 3.2

  params.my$tracks.inside <- 10
  params.my$tracks.outside <- 1

  params.my$line.width <- 1
  params.my$link.line.width <- 0.5

  params.my$text.size <-  0.6

  params.my$text.color <- 'black'

  params.my$track.padding <- c(0.07,  0.0, 0.07, 0.0,0.07, 0)

  params.my$grid.line.color <- 'lightgrey'
  params.my$chr.text.color <- 'grey'

  params.my$track.heights <- c(0.85, 0.07, 0.07, 0.1, 0.1,  0.1)
  params.my$track.height <- 0.1
  params.my$sub.tracks <- 1
  params.my$heatmap.cols <- c(alpha('lightcoral', 1),
                              alpha('lightcoral', 0.5),
                              alpha('lightgrey',0.10),
                              alpha('olivedrab2', 0.3),
                              alpha('olivedrab2', 0.5),
                              alpha('olivedrab2',.7),
                              alpha('olivedrab2', 0.75),
                              alpha('olivedrab3', 0.9),
                              alpha('olivedrab4', 0.9))
  params.my$heatmap.ranges <- c(0,1,3,4,8,16, 32,64,1000)

  #Set copynumber (and indel) colour scheme
  if (colour.scheme == 'picnic') {

    #tumour totalCN column
    params.my$heatmap.data.col.gain <- 6
    params.my$heatmap.data.col.loh <- 6

    params.my$heatmap.color.gain <- c(alpha('white',1), alpha('lightgrey',0.10), alpha('firebrick1',0.7), alpha('firebrick3',1.0))
    params.my$heatmap.ranges.gain <- c(0, 2, 4, 8, 1000)

    params.my$heatmap.ranges.loh <- c(0,1,2,1000)
    params.my$heatmap.color.loh <- c(alpha('darkblue', 1), alpha('grey',0.50), alpha('white',1))

    params.my$heatmap.key.gain.col <- alpha('firebrick1', 0.9)
    params.my$heatmap.key.loh.col <- alpha('darkblue', 1)
    params.my$heatmap.key.gain.title <- 'gain'
    params.my$heatmap.key.loh.title <- 'deletion'

    #Indel colours (to deferentiate from the copynumber colours)
    params.my$indel.mhomology <- 'mediumpurple4'
    params.my$indel.repeatmediated <- 'mediumpurple1'
    params.my$indel.other <- 'mediumorchid3'
    params.my$indel.insertion <- 'darkgreen'
    params.my$indel.complex <- 'grey'

  } else {
    #ascat
    params.my$heatmap.color.gain <- c( alpha('lightgrey',0.10), alpha('olivedrab2', 0.3),  alpha('olivedrab2', 0.5), alpha('olivedrab2',.7), alpha('olivedrab2', 0.75), alpha('olivedrab3', 0.9), alpha('olivedrab4', 0.9))
    params.my$heatmap.ranges.gain <- c(0,2,4,8,16, 32,64,1000)

    params.my$heatmap.ranges.loh <- c(0,1,1000)
    params.my$heatmap.color.loh <- c(alpha('lightcoral', 1), alpha('lightgrey',0.10))

    params.my$heatmap.key.gain.col <- alpha('olivedrab2', 0.3)
    params.my$heatmap.key.loh.col <- alpha('lightcoral', 1)
    params.my$heatmap.key.gain.title <- 'gain'
    params.my$heatmap.key.loh.title <- 'LOH'

    #tumour majorCN
    params.my$heatmap.data.col.gain <- 8
    #tumour minorCN
    params.my$heatmap.data.col.loh <- 7

    #Indel colours
    params.my$indel.mhomology <- 'firebrick4'
    params.my$indel.repeatmediated <- 'firebrick1'
    params.my$indel.other <- 'firebrick3'
    params.my$indel.insertion <- 'darkgreen'
    params.my$indel.complex <- 'grey'
  }

  return(params.my)

}

# main plotting function

#' Genome Plot
#'
#' Generates a plot for the visualisation of somatic variants across the genome, organised in a circle.
#' Variants plotted are single nucleotide variations (SNV), small insertions and deletions (indels),
#' copy number variations (CNV) and rearrangements.
#'
#' @param subsVcf.file SNV VCF file. The file should only contain SNV and should already be filtered according to the user preference, as all SNV in the file will be used and no filter will be applied.
#' @param indelsVcf.file Indels VCF file to be used to classify Indels and compute the proportion of indels at micro-homology. The files should only contain indels (no SNV) and should already be filtered according to the user preference, as all indels in the file will be used and no filter will be applied.
#' @param cnvsTab.file CNV TAB file (similar to ASCAT format). The file should be tab separated and contain a header in the first line with the following columns: 'seg_no', 'Chromosome', 'chromStart', 'chromEnd', 'total.copy.number.inNormal', 'minor.copy.number.inNormal', 'total.copy.number.inTumour', 'minor.copy.number.inTumour'
#' @param rearrBedpe.file SV (Rearrangements) BEDPE file. The file should contain a rearrangement for each row (two breakpoint positions should be on one row as determined by a pair of mates of paired-end sequencing) and should already be filtered according to the user preference, as all rearrangements in the file will be used and no filter will be applied. The file should contain a header in the first line with the following columns: "chrom1", "start1", "end1", "chrom2", "start2", "end2" and "sample" (sample name). In addition, either two columns indicating the strands of the mates, "strand1" (+ or -) and "strand2" (+ or -), or one column indicating the structural variant class, "svclass": translocation, inversion, deletion, tandem-duplication. The column "svclass" should correspond to (Sanger BRASS convention): inversion (strands +/- or -/+ and mates on the same chromosome), deletion (strands +/+ and mates on the same chromosome), tandem-duplication (strands -/- and mates on the same chromosome), translocation (mates are on different chromosomes).
#' @param sampleID Name of the sample.
#' @param genome.v set genome version: hg19, hg38, mm10 or canFam3.
#' @param file.ideogram name of the file that contain a user defined genome ideogram. Leave to NULL to load appropriate ideogram according to genome version.
#' @param plot_title title of the plot.
#' @param no_copynumber set to TRUE to disable plotting of copy number data
#' @param no_rearrangements set to TRUE to disable plotting of rearrangement data
#' @param no_indels set to TRUE to disable plotting of indels data
#' @param no_subs_legend set to TRUE to disable plotting of substitutions legend
#' @param out_format set to either png or svg
#' @param out_path directory where to place the plot.
#' @param rearr_only_assembled only include the rearrangements that have an assembly score
#' @param base.per.unit set RCircos base.per.unit parameter. Useful for whole exome data.
#' @return return the generated plot file name.
#' @export
genomePlot <- function(subsVcf.file, indelsVcf.file, cnvsTab.file, rearrBedpe.file,
                       sampleID, genome.v="hg19", ..., file.ideogram = NULL, plot_title = NULL,
                       no_copynumber = FALSE, no_rearrangements = FALSE, no_indels = FALSE,
                       no_subs_legend = FALSE, out_format = "png", out_path = ".",
                       rearr_only_assembled = FALSE, base.per.unit = NULL) {

  library(RCircos);
  library(scales)

  genome.bsgenome = switch(genome.v,
    "hg19" = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5,
    "hg38" = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
    "mm10" = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10,
    "canFam3" = BSgenome.Cfamiliaris.UCSC.canFam3::BSgenome.Cfamiliaris.UCSC.canFam3
    #"rn4"  = "BSgenome.Rnorvegicus.UCSC.rn4"
  )

  genome.ideogram = switch(genome.v,
    "hg19" = "UCSC.HG19.Human.CytoBandIdeogram",
    "hg38" = "UCSC.HG38.Human.CytoBandIdeogram",
    "mm10" = "UCSC.Mouse.GRCm38.CytoBandIdeogram",
    "canFam3" = "NA"
    #"rn4"  = "UCSC.Baylor.3.4.Rat.cytoBandIdeogram"
  )

  # load required genome and ideogram
  #library(genome.bsgenome, character.only = TRUE)

  if (!is.null(file.ideogram) && file.exists(file.ideogram)) {
    species.cyto <- read.table(file.ideogram, header=TRUE, sep='\t');
  } else {
    data(list=genome.ideogram,package = "RCircos");
    species.cyto <- get(genome.ideogram);
  }

  params.my <- set.plot.params()

  # rearrangement links colors
  inv.col <- alpha('dodgerblue2', 1)
  del.col <- alpha('coral2', 1) # originally .5
  #dupl.col <-  alpha('olivedrab3', 1) # originally 0.75
  dupl.col <-  alpha('darkgreen', 1)
  transloc.colour <- alpha('gray35', 1) # originally 0.5

  #Set up height, width and resolution parameters
  cPanelWidth = 0 #0.17
  graph.height = 4100
  graph.wd_ht_ratio = (5400/4100)  #width/height ratio
  graph.width = graph.height * graph.wd_ht_ratio
  graph.wd_res_ratio = (5400/550)
  graph.res = graph.width/graph.wd_res_ratio

  graph.height.inches = graph.height/graph.res
  graph.width.inches = graph.width/graph.res

  # indels
  indels <- NULL
  dels.formatted <- data.frame()
  ins.formatted <- data.frame()

  res <- vcfToIndelsClassification(indelsVcf.file, sampleID, genome.v)

  if (!no_indels && !is.null(res)) {
      indels <- res$indels_classified
      ins <- indels[which(indels$indel.type=='I'),]
      dels <- indels[which(indels$indel.type=='D' | indels$indel.type=='DI'),]
      ins$end <- ins$pos + ins$indel.length
      dels$end <- dels$pos + dels$indel.length
      ins.formatted <- ins[,c('chr', 'pos', 'end')];
      names(ins.formatted) <- c('Chromosome','chromStart','chromEnd')
      dels.formatted <- dels[,c('chr', 'pos', 'end')];
      names(dels.formatted) <- c('Chromosome','chromStart','chromEnd')
      if (nrow(ins.formatted)>0 && (genome.v=="hg19" || genome.v=="mm10")) {
        ins.formatted$Chromosome <- paste('chr',ins.formatted$Chromosome ,sep='')
      }
      if (nrow(dels.formatted)>0 && (genome.v=="hg19"|| genome.v=="mm10")) {
        dels.formatted$Chromosome <- paste('chr',dels.formatted$Chromosome ,sep='')
      }
      tile.cols <- vector()
      tile.cols[dels$classification=='Microhomology-mediated'] <- params.my$indel.mhomology #'firebrick4'
      tile.cols[dels$classification=='Repeat-mediated'] <- params.my$indel.repeatmediated #'firebrick1'
      tile.cols[dels$classification=='None'] <- params.my$indel.other #'firebrick3'
  }

  cat(paste(dim(indels)[1], ' indels \n'))

  # rearrangements
  rearrs.formatted <- data.frame()
  if (!no_rearrangements) {
    rearrs.formatted <- read.brass.bedpe(rearrBedpe.file, onlyAssembled = rearr_only_assembled)
    if(is.null(rearrs.formatted) || nrow(rearrs.formatted)==0){
      no_rearrangements <- TRUE
    }else{
      if (genome.v=="hg19" || genome.v=="mm10") {
        rearrs.formatted$Chromosome <- paste('chr', rearrs.formatted$Chromosome,sep='')
        rearrs.formatted$Chromosome.1 <- paste('chr', rearrs.formatted$Chromosome.1,sep='')
      }
    }
  }

  mut.order <- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T",
                 "G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T",
                 "A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T",
                 "G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T",
                 "A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T",
                 "G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T",
                 "A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T",
                 "G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T",
                 "A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T",
                 "G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T",
                 "A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T",
                 "G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")

  # substitutions
  subs.data <- getMutTables(subsVcf.file, onlyPASSED=FALSE, genome.v=genome.v, genomeSeq=genome.bsgenome,mut.order=mut.order)
  #subs.data <- getMutTablesTab(subsTab.file, onlyPASSED=FALSE, genomeSeq=get(genome.bsgenome))

  subs <- data.frame(chr=subs.data$muts$chroms,
                     position = subs.data$muts$starts,
                     wt = subs.data$muts$wt,
                     mt =  subs.data$muts$mt,
                     ref_base_pyrimidine_context =  subs.data$muts$pyrwt,
                     mutant_base_pyrimidine_context = subs.data$muts$pyrmut
                     )
  #subs <- subset(subs, chr!='Y')

  subs <- processSubs(subs)
  scatter.data <- subs$scatter.data
  scatter.colors <- subs$scatter.colors
  scatter.data.formatted <- data.frame(chromosome=scatter.data$chr,
                                       start=scatter.data$position,
                                       stop=scatter.data$position,
                                       logDistPrev=log10(scatter.data$distPrev),
                                       color=scatter.colors,
                                       mutType=scatter.data$mutType)
  cat(paste( dim(scatter.data)[1], ' subs \n'))

  # copy number
  cv.data <- data.frame()
  #Skip if no copynumber was requested
  if (!no_copynumber) {
    cv.data <- read.ascat(cnvsTab.file)
    if(is.null(cv.data) || nrow(cv.data)==0){
      no_copynumber <- TRUE
    } else if (genome.v=="hg19" || genome.v=="mm10" || genome.v=="canFam3") {
      cv.data$Chromosome <- paste('chr', cv.data$Chromosome,sep='')
    }
  }

################################################################################

  fn = file.path(out_path, paste(sampleID, ".genomePlot.", out_format, sep=''), fsep = .Platform$file.sep)

  if (out_format == 'png') {
    png(file=fn, height=graph.height, width=(graph.width*(1/(1-cPanelWidth))), res=graph.res)
  } else if (out_format == 'svg') {
    svg(fn, height=graph.height.inches, width=graph.width.inches)
  } else {
    stop("Invalid file type. Only png and svg are supported");
  }

  RCircos.Set.Core.Components(cyto.info=species.cyto, chr.exclude=NULL,  tracks.inside=params.my$tracks.inside, tracks.outside=params.my$tracks.outside);

  # set plot colours and parameters
  params <- RCircos.Get.Plot.Parameters();
  #RCircos doesn't allow resetting certain parameters
  #see implementation of RCircos.Reset.Plot.Parameters
  #params[which(names(params) %in% c("radius.len","plot.radius","chr.ideo.pos"))] <- NULL

  #toReset<-setdiff(names(params.my),c("radius.len","plot.radius","chr.ideo.pos","highlight.pos","chr.name.pos"))
  #params[toReset] <- params.my[toReset]

  #params$sub.tracks <- 1
  params$point.type <- params.my$point.type
  params$point.size <- params.my$point.size
  RCircos.Reset.Plot.Parameters(params)

  par(mar=c(0.001, 0.001, 0.001, 0.001))
  par(fig=c(cPanelWidth,0.75*(1-cPanelWidth)+cPanelWidth,0,1),cex=1.2)

  #RCircos.List.Parameters();

  #Set base.per.unit (especially needed when using exome file.ideogram)
  if (!is.null(base.per.unit)) {
    params <- RCircos.Get.Plot.Parameters();
    params$base.per.unit <- base.per.unit;
    RCircos.Reset.Plot.Parameters(params);
  }

  #RCircos.Set.Plot.Area(margins = 0.05);
  par(mai=c(params.my$genomeplot.margin, params.my$genomeplot.margin, params.my$genomeplot.margin, params.my$genomeplot.margin))
  plot.new()
  plot.window(c(-params.my$plot.radius,params.my$plot.radius), c(-params.my$plot.radius, params.my$plot.radius))
  RCircos.Chromosome.Ideogram.Plot.my(params.my$chr.text.color, params.my$grid.line.color, params.my$text.size);

  title(main = sampleID)

  if (!is.null(plot_title)) {
    title(paste(plot_title, sep=''), line=-1);
  }

  # substitutions
  if (exists("scatter.data.formatted")) {
    RCircos.Scatter.Plot.color(scatter.data=scatter.data.formatted , data.col=4, track.num=1, side="in", by.fold=0, scatter.colors = scatter.colors);
    cat('subs plotted \n')
  }

  #Set base.per.unit (especially needed when using exome file.ideogram)
  if (!is.null(base.per.unit)) {
    params <- RCircos.Get.Plot.Parameters();
    params$base.per.unit <- base.per.unit;
    RCircos.Reset.Plot.Parameters(params);
  }

  params <- RCircos.Get.Plot.Parameters();
  params$line.color <- 'white'
  params$highlight.width <- 0.2
  params$max.layers <- 5
  params$tile.color <- 'darkgreen'
  RCircos.Reset.Plot.Parameters(params)

  if (exists("ins.formatted") && nrow(ins.formatted)>0) {
      my.RCircos.Tile.Plot(tile.data=ins.formatted, track.num=5, side="in");
  }

  params <- RCircos.Get.Plot.Parameters();
  params$tile.color <- 'firebrick4'
  RCircos.Reset.Plot.Parameters(params)
  if (exists("dels.formatted") && nrow(dels.formatted)>0) {
      my.RCircos.Tile.Plot(tile.data=dels.formatted, track.num=6, side="in", tile.colors=tile.cols);
  }
  cat('indels plotted \n')

  # copy number
  if (exists('cv.data') && (nrow(cv.data)>0)) {

    heatmap.ranges.major <-params.my$heatmap.ranges.gain
    heatmap.color.major <-params.my$heatmap.color.gain
    heatmap.data.col.major <-params.my$heatmap.data.col.gain
    #heatmap.ranges.major <- c(0,2,4,8,16, 32,64,1000)
    #heatmap.color.major <- c( alpha('lightgrey',0.10), alpha('olivedrab2', 0.3),  alpha('olivedrab2', 0.5), alpha('olivedrab2',.7), alpha('olivedrab2', 0.75), alpha('olivedrab3', 0.9), alpha('olivedrab4', 0.9))

    RCircos.Heatmap.Plot.my(heatmap.data=cv.data, data.col=heatmap.data.col.major, track.num=7, side="in", heatmap.ranges=heatmap.ranges.major , heatmap.color=heatmap.color.major ); # major copy number

    #heatmap.ranges.minor <- c(0,1,1000)
    #heatmap.color.minor <- c(alpha('lightcoral', 1), alpha('lightgrey',0.10))
    heatmap.ranges.minor <-params.my$heatmap.ranges.loh
    heatmap.color.minor <-params.my$heatmap.color.loh
    heatmap.data.col.minor <-params.my$heatmap.data.col.loh
    RCircos.Heatmap.Plot.my(heatmap.data=cv.data, data.col=heatmap.data.col.minor, track.num=8, side="in", heatmap.ranges=heatmap.ranges.minor , heatmap.color=heatmap.color.minor ); # minor copy number

  }

  cat('copy number plotted \n')

  # rearrangements
  #Chromosome chromStart  chromEnd Chromosome.1 chromStart.1 chromEnd.1
  link.colors <- vector()
  if (exists("rearrs.formatted")) {
    link.data <- rearrs.formatted
    link.colors[link.data$pf==1 | link.data$pf==8] <- inv.col
    link.colors[link.data$pf==2] <- del.col
    link.colors[link.data$pf==4] <- dupl.col
    link.colors[link.data$pf==32] <- transloc.colour

    if (nrow( rearrs.formatted)>0) {
      RCircos.Link.Plot.my(link.data= rearrs.formatted, track.num=9,
                           by.chromosome=TRUE, link.colors);
      cat('rearrangements plotted \n')
    }
  }

  op <- par(lwd = 0.1)
  # square plotting region,
  # independent of device size
  ## At end of plotting, reset to previous settings:
  par(op)

  # side plots
  margins <- 0.25
  par(cex=0.5,mai = c(margins, margins, margins, margins))


  if (exists("subs.data") && !is.null(subs.data$passed.hist)) {
    # plot a histogram
    par(fig=c(cPanelWidth+0.74*(1-cPanelWidth),cPanelWidth+.995*(1-cPanelWidth),0.70, 0.95), new=TRUE)
    names(subs.data$passed.hist) <- NA
    barplot(subs.data$passed.hist , col=c(rep('royalblue',16), rep('black',16), rep('red', 16), rep('grey', 16), rep('green2', 16), rep('hotpink',16)),main=paste(nrow(subs.data$muts),'substitutions'), , border=NA)
  }

  if (!no_subs_legend) {
    par(fig=c(cPanelWidth+0.63*(1-cPanelWidth),cPanelWidth+0.72*(1-cPanelWidth),0.78,0.95), new=TRUE)
    # subs legend
    plot(1, type="n", axes=FALSE, xlab="", ylab="")
    legend(x='center',
         legend=c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'),
         col=c("royalblue", "black", "red", "grey", "green2", "hotpink"), pch=19, pt.cex=0.6,  horiz=FALSE, bty='n', )
  }

  par(fig=c(cPanelWidth+0.8*(1-cPanelWidth),cPanelWidth+.98*(1-cPanelWidth),0.46, 0.64), new=TRUE)

  # indels
  if (exists("indels") && !is.null(indels)) {
    indel.data <- c( sum(indels$indel.type=='D' & indels$classification == 'Microhomology-mediated'),
                      sum(indels$indel.type=='D' & indels$classification == 'Repeat-mediated'),
                      sum(indels$indel.type=='D' & indels$classification == 'None'),
                      sum(indels$indel.type=='I'),
                      sum(indels$indel.type=='DI'))
    #indel.col <- c('firebrick4', 'firebrick1', 'firebrick3', 'darkgreen', 'grey')
    indel.col <- c(params.my$indel.mhomology, params.my$indel.repeatmediated, params.my$indel.other, params.my$indel.insertion, params.my$indel.complex)
    indel.lbs <- c('deletion \n m-homology',  'deletion repeat', 'deletion other', 'insertion', 'complex')
    mp <- barplot(indel.data, main=paste(nrow(indels), 'deletions and insertions'), axes = FALSE, col=indel.col, axisnames = FALSE, width=1 , horiz=TRUE, border=NA)
    axis(2, at = mp, las=2,  labels = indel.lbs, col='grey', tick=FALSE, cex=0.5)
    axis(1, las=2, col='grey')
  }

  # copy number
  if (exists('cv.data') && (nrow(cv.data)>0)) {
    par(fig=c(cPanelWidth+0.73*(1-cPanelWidth),cPanelWidth+0.99*(1-cPanelWidth),0.32,0.42), new=TRUE) # copy number legend
    plot(1, type="n", axes=FALSE, xlab="", ylab="" , main=paste0('copy number'))
    legend(x='bottom',
         #legend=c('LOH', 'gain'),
         legend=c(params.my$heatmap.key.loh.title, params.my$heatmap.key.gain.title),
         #col=c(alpha('lightcoral', 1), alpha('olivedrab2', 0.3)), pch=15, pt.cex=2,  horiz=TRUE, bty='n',
         col=c(params.my$heatmap.key.loh.col, params.my$heatmap.key.gain.col), pch=15, pt.cex=2,  horiz=TRUE, bty='n',
         text.width=c(0.2,0.2))
  }

  if (exists("rearrs.formatted") && nrow( rearrs.formatted)) {
    par(fig=c(cPanelWidth+0.76*(1-cPanelWidth),cPanelWidth+ 0.98*(1-cPanelWidth),0.1,0.3), new=TRUE) # rearrangements

    rearrs.bar <- c(sum(rearrs.formatted$pf==32), sum(rearrs.formatted$pf==1)+sum(rearrs.formatted$pf==8), sum(rearrs.formatted$pf==2), sum(rearrs.formatted$pf==4))

    rearrs.col <- c(transloc.colour, inv.col, del.col, dupl.col)
    rearrs.lbs <- c('translocation',  'inversion', 'deletion', 't. duplication')
    mp <- barplot(rearrs.bar, main=paste(nrow(rearrs.formatted), 'rearrangements'), axes = FALSE, col=rearrs.col, axisnames = FALSE, width=1, horiz=TRUE, border=NA)
    axis(2, at = mp[rearrs.bar>0], las=2,  labels =rearrs.lbs[rearrs.bar>0], col='grey', tick=FALSE, cex=0.5)
    axis(1, las=2, col='grey')
  }

  dev.off();

  return(fn)
}

