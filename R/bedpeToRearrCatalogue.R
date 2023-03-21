
#' bedpeToRearrCatalogue
#'
#' This function converts a data frame BEDPE into a rearrangement catalogue,
#' you should pass rearrangements of only one sample, and one rearrangement for each paired-end mates.
#' The BEDPE data fram should contain the following columns: "chrom1", "start1", "end1", "chrom2", "start2", "end2" and "sample" (sample name).
#' In addition, either two columns indicating the strands of the mates, "strand1" (+ or -) and "strand2" (+ or -), or one column indicating the structural variant class, "svclass": translocation, inversion, deletion, tandem-duplication.
#' If you specify the "svclass" column, then the "strand1" and "strand2" columns will be ignored. If the "svclass" column is absent, then it will be created
#' using the convention of BrassII from the Sanger Institute pipeline: inversion when strand1 and strand2 are different, deletion when strand1 and strand2 are both +,
#' tandem-duplication when strand1 and strand2 are both -, and translocation when strand1 and strand2 are on different chromosomes.
#'
#' Please notice that the interpretation of strand1 and strand2 from your rearrangement caller may differ from the BrassII Sanger interpretation.
#' Typically, other callers may have the strand2 sign inverted with respect to what is shown by BrassII, so for example a deletion is when strand1 is + and strand2 is -,
#' instead of when both are + as in our case. To avoid confusion, double check the convention of your caller, and possibly specify the svclass column yourself to simply ignore
#' the strand1 and strand2 automated interpretation.
#' 
#' Optionally, the user can provide in the bedpe two additional columns, "non-template" and "micro-homology",
#' which should contain the DNA sequence inserted ("non-template") or deleted ("micro-homology") at the breakpoints junction.
#' A dot (".") should be inserted in these columns if a DNA sequence is not available. When these two columns are available,
#' a junctions catalogue will be computed and returned. A junction catalogue contains the counts of how many clustered/unclustered
#' rearrangements have non-templated insertions or micro-homology deletions of a certain size.
#'
#' @param sv_bedpe data frame BEDPE as described above
#' @param kmin minimum number of break points in a segment to consider it a cluster. Default is 10.
#' @param PEAK.FACTOR this factor is used to calculate a threshold for the minimum average distance of breakpoints in a cluster. The threshold is given by the expected distance divided by the PEAK.FACTOR. In turn, the expected distance is the number of base pairs in a genome divided by the total number of break points. Default is 10.
#' @return returns a list with the rearrangement catalogue (rearr_catalogue) for the given sample and the annotated bedpe (annotated_bedpe) for the given sample. Also, a junctions catalogue will be returned if the non-template and micro-homology columns are provided. If clusters of rearrangements are found then the clustering regions will also be returned (clustering_regions). 
#' @keywords bedpe rearrangement
#' @export
#' @examples
#' vcf_sv_file.bedpe <- "sample.bedpe"
#' sv_bedpe <- read.table(vcf_sv_file.bedpe,sep = "\t",header = TRUE,
#'                      stringsAsFactors = FALSE,check.names = FALSE)
#' #build a catalogue from the bedpe file
#' res <- bedpeToRearrCatalogue(sv_bedpe)
#' plotRearrSignatures(res$rearr_catalogue)
bedpeToRearrCatalogue <- function(sv_bedpe,
                                  kmin = 10,
                                  PEAK.FACTOR = 10){

  #check that the required columns are present
  required_cols <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2" , "sample")
  if(!length(intersect(required_cols,colnames(sv_bedpe)))==length(required_cols) & ("svclass" %in% colnames(sv_bedpe) | all(c("strand1","strand2") %in% colnames(sv_bedpe)))){
    stop("[error bedpeToRearrCatalogue] missing columns in subs data frame, following columns required: chrom1, start1, end1, chrom2, start2, end2 and sample. In addition either svclass or strand1 and strand2. Check ?bedpeToRearrCatalogue for details.")
  }

  clusters_table <- NULL
  
  if(nrow(sv_bedpe)>0){
    # make sure that there are no multiple samples here
    if(length(unique(sv_bedpe$sample))>1){
      # keep only the most frequent sample
      tmpsamplenames <- unique(sv_bedpe$sample)
      tmpcounts <- table(sv_bedpe$sample)
      pos <- which.max(tmpcounts)
      samplechoice <- names(pos)
      sv_bedpe <- sv_bedpe[sv_bedpe$sample==samplechoice,,drop=F]
      message("[warning bedpeToRearrCatalogue] bedpeToRearrCatalogue sample column should have only one sample name, however multiple sample names were detected: ",
              paste(tmpsamplenames,collapse = ", "),". Trying to fix this by using only the sample with the largest number of rearrangements (",samplechoice,"). ",
              "This fix should work if there are a few rearrangements from the germline sample, so all we would do is to remove the germline sample name and variants ",
              "while keeping all the tumour sample variants.")
    }
  }
  #Annotate the bedpe if necessary

  if(nrow(sv_bedpe)>0){

    clustering.result <- rearrangement.clustering_bedpe(sv_bedpe,
                                                        plot.path = NA,
                                                        kmin=kmin,
                                                        kmin.samples=1,
                                                        gamma.sdev=25,
                                                        PEAK.FACTOR=PEAK.FACTOR,
                                                        thresh.dist=NA)
    sv_bedpe <- clustering.result$sv_bedpe
    clustering_regions <- clustering.result$clustering_regions

    #check whether column svclass is present,
    #if not, compute it
    if (! "svclass" %in% colnames(sv_bedpe)){
      if ("strand1" %in% colnames(sv_bedpe) & "strand2" %in% colnames(sv_bedpe)){
        sv_bedpe <- classifyRearrangementsFromBedpe(sv_bedpe)
      }else{
        message("cannot classify rearrangements: svclass column missing, and cannot compute it because strand1 and strand2 are missing.")
      }
    }

    # remove translocations with length less than 1kb
    bkdist <- abs(sv_bedpe$start2 - sv_bedpe$start1)
    sv_bedpe[sv_bedpe$svclass!='translocation',"length"] <- bkdist[sv_bedpe$svclass!='translocation']
    all_sv_annotated <- sv_bedpe
    toberemoved <- sv_bedpe$svclass!='translocation' & bkdist<1e3
    all_sv_annotated[toberemoved,"FILTER"] <- "length<1e3"
    all_sv_annotated[!toberemoved,"FILTER"] <- "PASS"
    sv_bedpe <- sv_bedpe[!toberemoved,]


  }else{

    all_sv_annotated <- NULL
    clustering_regions <- NULL
  }
  
  #now compute the catalogue (this takes care of the 0 SVs case)
  rearr_catalogue <- prepare.rearr.catalogue_fromAnnotatedBedpe(sv_bedpe)
  # get the junctions catalogue too (this takes care of the 0 SVs case)
  junctions_catalogue <- build_junctions_catalogue(sv_bedpe)
  if(!is.null(junctions_catalogue)) colnames(junctions_catalogue) <- colnames(rearr_catalogue)

  # return all results
  return_list <- list()
  return_list$rearr_catalogue <- rearr_catalogue
  return_list$junctions_catalogue <- junctions_catalogue
  return_list$annotated_bedpe <- all_sv_annotated
  if(!is.null(clustering_regions)){
    if(nrow(clustering_regions)>0) return_list$clustering_regions <- clustering_regions
  }
  return(return_list)
}


# source('../lib/pcf/fastPCF.R')
# source('../lib/pcf/extract.kat.regions.R')
# source('../lib/utils/calcIntermutDist.R')

# this is used for per-sample clustering of both single-base substitutions and rearrangement breakpoints

rearrangement.clustering_bedpe <- function(sv_bedpe,
                                           plot.path=NA,
                                           kmin=10,# how many points at minimum in a peak, for the pcf algorithm
                                           kmin.samples=kmin, # how many different samples at minimum in  a peak
                                           gamma.sdev=25, #
                                           PEAK.FACTOR=4,
                                           thresh.dist=NA,
                                           gamma=NA,
                                           kmin.filter=kmin # if the pcf parameter is different from the definition of a peak
) {

  #add an id to the rearrangement
  sv_bedpe$id <- 1:nrow(sv_bedpe)

  #functions below expect rows to be organised by chromosomes and ordered by position on the chromosome

  #prepare a dataframe for the calculation
  rearrs.left <- sv_bedpe[,c('chrom1','start1','sample')]
  names(rearrs.left ) <- NA
  rearrs.right <- sv_bedpe[,c('chrom2','start2','sample')]
  names(rearrs.right ) <- NA
  rearrs.cncd <- rbind(rearrs.left , rearrs.right  )
  colnames(rearrs.cncd) <- c('chr', 'position', 'sample')
  rearrs.cncd$isLeft <- c(rep(TRUE, nrow(rearrs.left)), rep(FALSE, nrow(rearrs.left)))
  rearrs.cncd$id <- c(sv_bedpe$id, sv_bedpe$id)
  # sample.bps <- rearrs.cncd
  #need to reorder
  sample.bps <- NULL
  for (chrom_i in unique(rearrs.cncd$chr)){
    tmptab <- rearrs.cncd[rearrs.cncd$chr==chrom_i,,drop=FALSE]
    tmptab <- tmptab[order(tmptab$position),,drop=FALSE]
    sample.bps <- rbind(sample.bps,tmptab)
  }
  rownames(sample.bps) <- 1:nrow(sample.bps)

  #run the algorithm
  genome.size <- 3 * 10^9
  MIN.BPS <- 10 # minimal number of breakpoints on a chromosome to do any any segmentation

  logScale <- FALSE

  exp.dist <-genome.size/nrow(sample.bps)

  if (logScale) {
    sample.bps$intermut.dist <- log10(calcIntermutDist(sample.bps, first.chrom.na=FALSE)$distPrev) # calculate the distances between the breakpoints
    if (is.na(thresh.dist)) {
      thresh.dist <- log10(exp.dist/PEAK.FACTOR) # calculate the threshold to call a peak
    }
  } else {

    sample.bps$intermut.dist <- calcIntermutDist(sample.bps, first.chrom.na=FALSE)$distPrev
    if (is.na(thresh.dist)) {
      thresh.dist <- exp.dist/PEAK.FACTOR
    }
  }


  if (is.na(gamma) & !is.na(gamma.sdev)) {
    # compute the mean absolute deviation
    sdev <- getMad(sample.bps$intermut.dist);
    gamma <- gamma.sdev*sdev
  }



  sample.bps$is.clustered.single <- rep(FALSE, nrow(sample.bps))

  all.kat.regions <- data.frame()

  for (chrom in unique(sample.bps$chr)) { # loop over chromosomes

    sample.bps.flag <- sample.bps$chr==chrom #   breakpoints on a current chromosome
    # sample.bps.chrom <- sample.bps[sample.bps.flag,]
    # sample.bps.chrom <- sample.bps.chrom[order(sample.bps.chrom$position),]
    #
    if (sum(sample.bps.flag )>MIN.BPS ) { # if there are enough breakpoints on a chromosome to run pcf

      data.points <- sample.bps$intermut.dist[sample.bps.flag]
      # data.points <- sample.bps.chrom$intermut.dist

      res = exactPcf(data.points, kmin, gamma, T)

      #reorder results
      sample.bps$mean.intermut.dist[sample.bps.flag] <- res$yhat

      # prepare the points for pcf
      subs <- data.frame(chr=sample.bps$chr[sample.bps.flag], pos=sample.bps$position[sample.bps.flag], sample=sample.bps$sample[sample.bps.flag])
      kat.regions <- extract.kat.regions(res, thresh.dist, subs, doMerging=TRUE, kmin.samples=1,  kmin.filter= kmin.filter) # extract peaks, this is special case as we want at least kmin samples

      all.kat.regions <- rbind(all.kat.regions, kat.regions)
      if (!is.null(kat.regions) && nrow( kat.regions )>0) { # if there are any kataegis regions found on this chormosome
        for (k in 1:nrow(kat.regions)) {

          sample.bps$is.clustered.single[which(sample.bps.flag)[ kat.regions$firstBp[k] : kat.regions$lastBp[k]]] <- TRUE # call all breakpoints as clustered
        }
      }
    } else {

      sample.bps$mean.intermut.dist[sample.bps.flag] <- mean(sample.bps$intermut.dist[sample.bps.flag])
    }
  }



  if (!logScale) { # even if pcf was run on non-logged distances, I log the output
    sample.bps$intermut.dist <- log10(sample.bps$intermut.dist)
    sample.bps$mean.intermut.dist <- log10(sample.bps$mean.intermut.dist)
  }

  # a rearrangement is in a cluster if any of its breakpoints are
  sample.bps$is.clustered <- sample.bps$is.clustered.single
  sample.bps$is.clustered[sample.bps$id %in% subset(sample.bps, is.clustered.single==TRUE)$id] <- TRUE

  # mark both breakpoints of a rearrangement as clustered if any is
  sv_bedpe$is.clustered <- sv_bedpe$id %in% sample.bps$id[sample.bps$is.clustered]

  result <- list()
  result$sv_bedpe <- sv_bedpe
  result$clustering_regions <- all.kat.regions
  result
}


prepare.rearr.catalogue_fromAnnotatedBedpe <- function(sv_bedpe) {

  catalogue.labels <- c('clustered_del_1-10Kb', 'clustered_del_10-100Kb', 'clustered_del_100Kb-1Mb', 'clustered_del_1Mb-10Mb', 'clustered_del_>10Mb', 'clustered_tds_1-10Kb', 'clustered_tds_10-100Kb', 'clustered_tds_100Kb-1Mb', 'clustered_tds_1Mb-10Mb', 'clustered_tds_>10Mb', 'clustered_inv_1-10Kb', 'clustered_inv_10-100Kb', 'clustered_inv_100Kb-1Mb', 'clustered_inv_1Mb-10Mb', 'clustered_inv_>10Mb', 'clustered_trans', 'non-clustered_del_1-10Kb', 'non-clustered_del_10-100Kb', 'non-clustered_del_100Kb-1Mb', 'non-clustered_del_1Mb-10Mb', 'non-clustered_del_>10Mb', 'non-clustered_tds_1-10Kb', 'non-clustered_tds_10-100Kb', 'non-clustered_tds_100Kb-1Mb', 'non-clustered_tds_1Mb-10Mb', 'non-clustered_tds_>10Mb', 'non-clustered_inv_1-10Kb', 'non-clustered_inv_10-100Kb', 'non-clustered_inv_100Kb-1Mb', 'non-clustered_inv_1Mb-10Mb', 'non-clustered_inv_>10Mb', 'non-clustered_trans')

  all_catalogues <- as.data.frame(matrix(nrow = length(catalogue.labels),ncol = 0))
  rownames(all_catalogues) <- catalogue.labels

  if (nrow(sv_bedpe)>0){
    for (sample_name in unique(sv_bedpe$sample)){
      sample.rearrs <- sv_bedpe[sv_bedpe$sample==sample_name,]

      rearr_catalogue <- as.data.frame(matrix(0,nrow = length(catalogue.labels),ncol = 1))

      if (nrow(sample.rearrs)>0) {

        label1 <- rep('non-clustered', nrow(sample.rearrs))
        label1[ sample.rearrs$is.clustered] <- 'clustered'

        label2 <- rep('', nrow(sample.rearrs))
        label2[ sample.rearrs$svclass=='deletion'] <- '_del'
        label2[ sample.rearrs$svclass=='translocation'] <- '_trans'
        label2[ sample.rearrs$svclass=='inversion'] <- '_inv'
        label2[ sample.rearrs$svclass=='tandem-duplication'] <- '_tds'

        label3 <- rep('', nrow(sample.rearrs))
        sample.rearrs$bkdist <- abs(sample.rearrs$start2 - sample.rearrs$start1)
        label3[ sample.rearrs$svclass!='translocation' & sample.rearrs$bkdist<=1e4] <- '_1-10Kb'
        label3[ sample.rearrs$svclass!='translocation' & sample.rearrs$bkdist>1e4 & sample.rearrs$bkdist<=1e5 ] <- '_10-100Kb'
        label3[ sample.rearrs$svclass!='translocation' & sample.rearrs$bkdist>1e5 & sample.rearrs$bkdist<=1e6 ] <- '_100Kb-1Mb'
        label3[ sample.rearrs$svclass!='translocation' & sample.rearrs$bkdist>1e6 & sample.rearrs$bkdist<=1e7 ] <- '_1Mb-10Mb'
        label3[ sample.rearrs$svclass!='translocation' & sample.rearrs$bkdist>1e7 ] <- '_>10Mb'

        sample.rearrs$catalogue.label <- paste0(label1, label2, label3)

        sample.table <- as.data.frame(table( sample.rearrs$catalogue.label),drop=FALSE)
        rownames(sample.table ) <- sample.table$Var

        rearr_catalogue <-  sample.table [as.character(catalogue.labels), 'Freq',drop=FALSE ]

      }

      rearr.catalogue <- rearr_catalogue
      rownames(rearr.catalogue) <- catalogue.labels
      colnames(rearr.catalogue) <- sample_name
      rearr.catalogue[is.na(rearr.catalogue)] <- 0

      all_catalogues <- cbind(all_catalogues,rearr.catalogue)
    }
  }else{
    all_catalogues <- as.data.frame(matrix(0,nrow = length(catalogue.labels),ncol = 1))
    rownames(all_catalogues) <- catalogue.labels
  }
  all_catalogues
}



classifyRearrangementsFromBedpe <- function(sv_bedpe){
  svclass <- c()
  for (i in 1:nrow(sv_bedpe)){
    if(sv_bedpe[i,"chrom1"]!=sv_bedpe[i,"chrom2"]){
      svclass <- c(svclass,"translocation")
    }else if(sv_bedpe[i,"strand1"]!=sv_bedpe[i,"strand2"]){
      svclass <- c(svclass,"inversion")
    }else if(sv_bedpe[i,"strand1"]=="+"){
      svclass <- c(svclass,"deletion")
    }else if(sv_bedpe[i,"strand1"]=="-"){
      svclass <- c(svclass,"tandem-duplication")
    }
  }
  sv_bedpe[,"svclass"] <- svclass
  #return updated df
  sv_bedpe
}


build_junctions_catalogue <- function(annotated_bedpe){
  
  junctions_catalogue_channels <- paste(rep(c("clustered","non-clustered"),each=7),
                                        rep(c(rep("_non-templated",3),rep("_homologous",3),"_other"),2),
                                        rep(c(rep(c("_1-3","_4-10","_>10"),2),""),2),sep = "")
  junctions_catalogue <- data.frame(sample=rep(0,length(junctions_catalogue_channels)),
                                    row.names = junctions_catalogue_channels,
                                    stringsAsFactors = F)
  if(nrow(annotated_bedpe)==0){
    # cut short here returning a 0 catalogue if there are no SVs
    return(junctions_catalogue)
  }
  
  # # check if non-template and micro-homology columns are present
  if(all(c("non-template","micro-homology") %in% colnames(annotated_bedpe))){
    for (clustered in c(TRUE,FALSE)){
      # clustered <- F
      channelc <- ifelse(clustered,"clustered","non-clustered")
      for (typebp in c("_non-templated","_homologous","_other")){
        # typebp <- "_non-templated"
        channel <- paste0(channelc,typebp)
        if(typebp %in% c("_non-templated","_homologous")){
          bedpecol <- ifelse(typebp=="_non-templated","non-template","micro-homology")
          currentseqlength <- nchar(annotated_bedpe[annotated_bedpe[,bedpecol]!="." & annotated_bedpe$is.clustered==clustered,bedpecol])
          if(length(currentseqlength)>0){
            currentlengthtable <- table(currentseqlength)
            for(i in 1:length(currentlengthtable)){
              # i <- 1
              tmpn <- as.numeric(names(currentlengthtable))[i]
              if(tmpn>0 & tmpn<=3){
                junctions_catalogue[paste0(channel,"_1-3"),"sample"] <- junctions_catalogue[paste0(channel,"_1-3"),"sample"] + currentlengthtable[i]
              }else if(tmpn>3 & tmpn<=10){
                junctions_catalogue[paste0(channel,"_4-10"),"sample"] <- junctions_catalogue[paste0(channel,"_4-10"),"sample"] + currentlengthtable[i]
              }else if(tmpn>10){
                junctions_catalogue[paste0(channel,"_>10"),"sample"] <- junctions_catalogue[paste0(channel,"_>10"),"sample"] + currentlengthtable[i]
              }
            }
          }
        }else{
          tmpn <- sum(annotated_bedpe[,"non-template"]=="." &  annotated_bedpe[,"micro-homology"]=="." & annotated_bedpe$is.clustered==clustered)
          junctions_catalogue[channel,"sample"] <- junctions_catalogue[channel,"sample"] + tmpn
        }
      }
    }
    return(junctions_catalogue)
  }else{
    # return NULL if there are SVs but no "non-template" and "micro-homology" columns
    return(NULL)
  }
  
}
