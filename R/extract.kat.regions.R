# source('../lib/utils/assignPvalues.R')

extract.kat.regions <- function (res, imd, subs,  kmin.samples=10, pvalue.thresh=1, rate.factor.thresh=1, doMerging=FALSE, kmin.filter=NA, bp.rate=NA) {
	
	segInterDist <-  res$yhat
	kataegis.threshold <- imd
	
	kat.regions.all = data.frame()	
	
	chr <- as.character(subs$chr[1])
        
	positions <- subs$pos

	katLoci = (segInterDist<=kataegis.threshold) # flag specifying if a point is in a peak
        
        if(sum(katLoci)>0) {
            
            start.regions = which(katLoci[-1] & !(katLoci[-(length(katLoci))]) # katLoci breakpoints			
                | (katLoci[-1] & katLoci[-(length(katLoci))] & segInterDist[-1] != segInterDist[-length(katLoci)] )
                                  )+1 # endpoints between peaks
            if (katLoci[1]) {start.regions <- c(1, start.regions)}
            
            end.regions = which(!(katLoci[-1]) & katLoci[-(length(katLoci))] #
                | (katLoci[-1] & katLoci[-(length(katLoci))] & segInterDist[-1] != segInterDist[-length(katLoci)] )
                                ) #
            if (katLoci[length(katLoci)]) {end.regions <- c( end.regions, length(katLoci))}
            
            start.regions.init <- start.regions
            end.regions.init <- end.regions 
            
                                        # handling special cases
            if(length(end.regions)+length(start.regions)>0) {  # if there are any discontinuities in the segmentation at all 						
                if (length(end.regions)==1 & length(start.regions)==0){ 
                    start.regions <- 1                                    
                } else if (length(start.regions)==1 & length(end.regions)==0){                                    
                    end.regions <- length(positions)                                    
                } else if ((end.regions[1]<start.regions[1])&& (start.regions[length(start.regions)]>end.regions[length(end.regions)])) {
                                        # starts and ends are the same length, but missing both endpoints
                                    
                    start.regions <- c(1,start.regions)
                    end.regions <- c(end.regions,  length(positions))
                    
                } else if (end.regions[1]<start.regions[1]){
                                        # starts will be one shorter
                    start.regions <- c(1, start.regions)
                    
                } else if (start.regions[length(start.regions)]>end.regions[length(end.regions)]){
                                        # ends will be one shorter
                    
                    end.regions <- c(end.regions,  length(positions))
                }
                
                if (length(start.regions)!=length(end.regions)) {
                    browser()
                }

                               

                                # prepare a data structure that will be later filled up
                kat.regions.all <- data.frame(
                    chr=subs$chr[1],
                    start.bp=rep(NA,length(start.regions)), # start coordinate [bp]
                    end.bp=rep(NA,length(start.regions)), # end coordinate [bp]
                    length.bp=rep(NA,length(start.regions)), # length [bp]
                    number.bps=rep(NA,length(start.regions)),
                    number.bps.clustered=rep(NA,length(start.regions)),
                    avgDist.bp=rep(NA,length(start.regions)),
                    no.samples=rep(NA,length(start.regions)),
                    no.del =rep(NA,length(start.regions)),
                    no.dup =rep(NA,length(start.regions)),
                    no.inv= rep(NA,length(start.regions)),
                    no.trn = rep(NA,length(start.regions)),
                    firstBp=start.regions,
                    lastBp=end.regions                                    )
                
                kat.regions.all <- hotspotInfo(kat.regions.all, subs, segInterDist)

                step.segInterDist.left <- rep(NA, length(segInterDist))
                step.segInterDist.left[2:length(segInterDist)] <- segInterDist[2:length(segInterDist)]- segInterDist[1:(length(segInterDist)-1)]       
                step.segInterDist.right <- rep(NA, length(segInterDist))
                step.segInterDist.right[1:(length(segInterDist)-1)] <- segInterDist[1:(length(segInterDist)-1)]- segInterDist[2:(length(segInterDist))]

                kat.regions.all$step.left <-  step.segInterDist.left[start.regions]
                kat.regions.all$step.right <-  step.segInterDist.right[end.regions]
               

                                        # run the filters on the regions of increased frequency
                                        # make sure there are at least kmin samples
                
                if ((!is.null(kat.regions.all)) && (nrow(kat.regions.all)>0)) {
                    kat.regions.all <- subset(kat.regions.all, no.samples>=kmin.samples)
                }
                
         
                                # make sure there are at least kmin.filter breakpoints
                if (!is.na(kmin.filter)) {
                    kat.regions.all <- subset(kat.regions.all, number.bps>=kmin.filter)
                }

                
                                
                                        # make sure the p-value is less than somethng
                if ((!is.null(kat.regions.all)) && (nrow(kat.regions.all)>0)) {
                    kat.regions.all <- assignPvalues(kat.regions.all, subs, bp.rate=bp.rate)
                    kat.regions.all <- subset(kat.regions.all, pvalue<=pvalue.thresh)
                                        # only keep the hotspots that exceed the theshold
                    kat.regions.all <- subset(kat.regions.all, rate.factor>=rate.factor.thresh)
                }  
                
                                       
                                
                                
                                
                                        # merge segments if both were found to be peaks
                if (doMerging) {
                    if(nrow(kat.regions.all)>1){
                        for(r in 2:nrow(kat.regions.all)){
                            if (kat.regions.all$lastBp[r-1] == (kat.regions.all$firstBp[r]-1)) {
                                        # merge two segments
                                kat.regions.all$firstBp[r] <- kat.regions.all$firstBp[r-1]
                                kat.regions.all$firstBp[r-1] <- NA
                                kat.regions.all$lastBp[r-1] <- NA
                                kat.regions.all$avgDist.bp[r] <- NA # this will need to be updated as segments are being merged
                            }
                        }        
                    }
                                        # remove some of the merged segments
                    kat.regions.all <- subset(kat.regions.all, !is.na(firstBp) & !is.na(lastBp))
                    
                                        # update the info on hotspots that might have changed when they were merged
                    kat.regions.all <- hotspotInfo( kat.regions.all ,  subs, segInterDist)
                    kat.regions.all <- assignPvalues(kat.regions.all, subs, bp.rate=bp.rate)
                } # end merging
                
                                
                
                
            } # end if there are discontinuities in the segmentation
        } # if there are any points under the inter-mutation distance threshold
        
	kat.regions.all
        
    }

hotspotInfo <- function(kat.regions.all, subs, segInterDist=c()) {
    if(nrow(kat.regions.all)>0){
        for(r in 1:nrow(kat.regions.all)){
            
                                        # indices of the breakpoints in the hotspot
            subs.hotspot <-subs[kat.regions.all$firstBp[r]:kat.regions.all$lastBp[r],]
            
            kat.regions.all[r,'start.bp'] <- min(subs.hotspot$pos)
            kat.regions.all[r,'end.bp'] <- max(subs.hotspot$pos)
            kat.regions.all[r,'length.bp'] <-  kat.regions.all[r,'end.bp'] - kat.regions.all[r,'start.bp'] 
            kat.regions.all[r,'number.bps'] <- nrow(subs.hotspot)
            kat.regions.all[r,'number.bps.clustered'] <- sum(subs.hotspot$is.clustered)
            
            if (length(segInterDist)>0 & is.na(kat.regions.all[r,'avgDist.bp'])) {
                kat.regions.all[r,'avgDist.bp'] <- mean(segInterDist[kat.regions.all$firstBp[r]:kat.regions.all$lastBp[r]])
            }
            kat.regions.all[r,'no.samples'] <- length(unique(subs.hotspot$sample))

            if ('pf' %in% colnames(subs.hotspot)){
                kat.regions.all[r,'no.del'] <- nrow(subset(subs.hotspot, pf==2))
                kat.regions.all[r,'no.dup'] <- nrow(subset(subs.hotspot, pf==4))
                kat.regions.all[r,'no.inv'] <- nrow(subset(subs.hotspot, pf==1 | pf==8))
                kat.regions.all[r,'no.trn'] <- nrow(subset(subs.hotspot, pf==32))
            }
            
        } # for all peaks
    } # if there is at least one peak
    kat.regions.all
}
