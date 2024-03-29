#' find kataegis
#' 
#' Given a table with a list of single nucleotide variants (SNVs) this function
#' identifies groups of SNVs that are narrowly clustered together, and that are
#' likely part of a kataegis. An annotated table of SNVs and a list of kataegis
#' regions are returned.
#' 
#' @param snvs_table table listing the single nucleotide variants to use. Necessary column names are: chr, position, REF, ALT.
#' @param sample_name name of sample
#' @return kataegis regions and annotated SNVs
#' @export
findKataegis <- function(snvs_table,
                         sample_name){
  
  subs.bps <- data.frame(chr=as.character(snvs_table$chr ),
                         position=snvs_table$position,
                         pf=-1,
                         sample=sample_name,
                         id=as.character(1:nrow(snvs_table)),
                         stringsAsFactors = F)
  
  subs_clustering_result <- subs_clustering(subs.bps,
                                            kmin=2,# two mutations
                                            kmin.samples=1, # they can belong to one sample only
                                            gamma.sdev=25,
                                            PEAK.FACTOR=NA,
                                            thresh.dist=1000,
                                            kmin.filter=6)
  
  
  # save the results
  result <- list()
  result$katregions <- subs_clustering_result$kat.regions
  snvs_table$is.kataegis <- subs_clustering_result$sample.bps$is.clustered
  result$snvs_table <- snvs_table
  return(result)
}

# this is used for per-sample clustering of both single-base substitutions and rearrangement breakpoints

subs_clustering <- function(sample.bps,
                            kmin=10,# how many points at minimum in a peak, for the pcf algorithm
                            kmin.samples=kmin, # how many different samples at minimum in  a peak
                            gamma.sdev=25, #
                            PEAK.FACTOR=4,
                            thresh.dist=NA,
                            gamma=NA,
                            kmin.filter=kmin) { # if the pcf parameter is different from the definition of a peak

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
    
    
    if (sum(sample.bps.flag )>MIN.BPS ) { # if there are enough breakpoints on a chromosome to run pcf
      
      data.points <- sample.bps$intermut.dist[sample.bps.flag]
      
      
      res = exactPcf(data.points, kmin, gamma, T)
      
      sample.bps$mean.intermut.dist[sample.bps.flag] <- res$yhat
      
      # prepare the points for pcf
      subs <- data.frame(
        chr=sample.bps$chr[sample.bps.flag],
        pos=sample.bps$position[sample.bps.flag],
        sample=sample.bps$sample[sample.bps.flag],
        id=sample.bps$id[sample.bps.flag]
      )
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
  
  result <- list()
  result$sample.bps <- sample.bps
  result$kat.regions <- all.kat.regions
  result
}
