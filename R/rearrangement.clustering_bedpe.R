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
  
  #prepare a dataframe for the calculation
  rearrs.left <- sv_bedpe[,c('chrom1','start1','sample','name')];
  names(rearrs.left ) <- NA
  rearrs.right <- sv_bedpe[,c('chrom2','start2','sample','name')];
  names(rearrs.right ) <- NA
  rearrs.cncd <- rbind(rearrs.left , rearrs.right  );
  rearrs.cncd$isLeft <- c(rep(TRUE, nrow(rearrs.left)), rep(FALSE, nrow(rearrs.left)))
  colnames(rearrs.cncd) <- c('chr', 'position', 'sample', 'id')
  sample.bps <- rearrs.cncd

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
      
      
      if (sum(sample.bps.flag )>MIN.BPS ) { # if there are enough breakpoints on a chromosome to run pcf
          
          data.points <- sample.bps$intermut.dist[sample.bps.flag]

          
          res = exactPcf(data.points, kmin, gamma, T)

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
  
  cat(paste('plotting to', plot.path))
  if (!is.na(plot.path)) {
      plotScatterCirco(sample.bps, plot.path)
  } 
  
  # mark both breakpoints of a rearrangement as clustered if any is
  sv_bedpe$is.clustered <- sv_bedpe$name %in% sample.bps$id[sample.bps$is.clustered.single]
  
  result <- list()
  result$sv_bedpe <- sv_bedpe
  result$kat.regions <- all.kat.regions
  result
}


