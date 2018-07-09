#part of the code for clustering rearrangements
assignPvalues <- function(kat.regions, chrom.bps, bp.rate=NA) {
  
  if (is.na(bp.rate)) { # estimate the chromosome rate
    left.bp <- min(chrom.bps$pos)
    right.bp <-  max(chrom.bps$pos)
    bp.rate <- nrow(chrom.bps)/ (right.bp - left.bp)
  }
  
  # assume binomial distribution
  kat.regions$pvalue <- 1-pbinom(kat.regions$number.bps, kat.regions$end.bp - kat.regions$start.bp, bp.rate)
  
  kat.regions$d.seg<- (kat.regions$number.bps/( kat.regions$end.bp - kat.regions$start.bp))
  
  kat.regions$rate.factor <- kat.regions$d.seg/bp.rate
  
  kat.regions
}

calcIntermutDist <- function (subs.type, first.chrom.na = FALSE) {
  
  subs.type.processed <- data.frame()
  for (c in unique(subs.type$chr)) {
    # choose subs from only one chromosome at a time
    
    subs.type.chrom <- subset(subs.type, subset=subs.type$chr==c)
    # sort the subs by position
    subs.type.chrom <- subs.type.chrom [order(subs.type.chrom$position),]
    
    if (first.chrom.na) {
      subs.type.chrom$prevPos <- c(NA,subs.type.chrom$position[1:nrow(subs.type.chrom)-1])
    } else {
      subs.type.chrom$prevPos <- c(0,subs.type.chrom$position[1:nrow(subs.type.chrom)-1])        
    }
    subs.type.chrom$distPrev  <- subs.type.chrom$position -  subs.type.chrom$prevPos
    
    subs.type.processed <- rbind(subs.type.processed,subs.type.chrom)
  }
  
  subs.type.processed$distPrev[subs.type.processed$distPrev==0] <- 1
  subs.type.processed 
}



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


# source("../lib/merge.with.order.R")

generateHist <- function(mutTable, normalise=TRUE, pyr.transcribed=NA,mut.order) {
  
  if(nrow(mutTable)>0) {
    
    if (! 'mut.context' %in% names(mutTable)) {
      mutTable$mut.context <- paste(mutTable$pyrbbef, '[',mutTable$pyrwt, '>',mutTable$pyrmut,']', mutTable$pyrbaft, sep='')
    }
    
    
    if (length(pyr.transcribed)==1) {
      mut.context.table <- as.data.frame(table(mutTable$mut.context))
      rownames(mut.context.table) <- as.character(mut.context.table$Var1)
      mut.context.table <- merge.with.order(data.frame(mut=mut.order),  mut.context.table , by.x='mut', by.y='Var1', all.x=TRUE, keep_order=TRUE)
      mut.context.table$Freq[is.na(mut.context.table$Freq)] <- 0
      total.muts <- sum(mut.context.table[,'Freq'])
      
      if (normalise) {
        signHist <- mut.context.table[,'Freq']/total.muts
      } else {
        signHist <- mut.context.table[,'Freq']
      }
      
      
    } else {
      
      mut.context.table.transcribed <- as.data.frame(table(mutTable$mut.context[pyr.transcribed]))
      rownames(mut.context.table.transcribed) <- as.character(mut.context.table.transcribed$Var1)
      mut.context.table.transcribed <- merge.with.order(data.frame(mut=mut.order), mut.context.table.transcribed , by.x='mut', by.y='Var1', all.x=TRUE, , keep_order=TRUE)
      mut.context.table.transcribed$Freq[is.na(mut.context.table.transcribed$Freq)] <- 0
      
      mut.context.table.nontranscribed <- as.data.frame(table(mutTable$mut.context[!pyr.transcribed]))
      rownames(mut.context.table.nontranscribed) <- as.character(mut.context.table.nontranscribed$Var1)
      mut.context.table.nontranscribed <- merge.with.order(data.frame(mut=mut.order), mut.context.table.nontranscribed , by.x='mut', by.y='Var1', all.x=TRUE, , keep_order=TRUE)
      mut.context.table.nontranscribed$Freq[is.na(mut.context.table.nontranscribed$Freq)] <- 0
      
      total.muts <- sum(mut.context.table.transcribed[,'Freq'])+ sum(mut.context.table.nontranscribed[,'Freq'])
      
      signHist <- rbind(mut.context.table.transcribed[,'Freq'], mut.context.table.nontranscribed[,'Freq'])
      if (normalise) {
        signHist <-sign.hist/total.muts
      }
      
    }
    
    
    
  } else {
    signHist <- rep(0,96)
    
  }
  
  names(signHist) <- mut.order
  
  
  signHist
}




merge.with.order <- function(x,y, ..., sort = T, keep_order)
{
  
  # this function works just like merge, only that it adds the option to return the merged data.frame ordered by x (1) or by y (2)
  add.id.column.to.data <- function(DATA)
  {
    data.frame(DATA, id... = seq_len(nrow(DATA)))
  }
  # add.id.column.to.data(data.frame(x = rnorm(5), x2 = rnorm(5)))
  order.by.id...and.remove.it <- function(DATA)
  {
    # gets in a data.frame with the "id..." column.  Orders by it and returns it
    if(!any(colnames(DATA)=="id...")) stop("The function order.by.id...and.remove.it only works with data.frame objects which includes the 'id...' order column")
    
    ss_r <- order(DATA$id...)
    ss_c <- colnames(DATA) != "id..."
    DATA[ss_r, ss_c]		
  }
  
  # tmp <- function(x) x==1; 1	# why we must check what to do if it is missing or not...
  # tmp()
  
  if(!missing(keep_order))
  {
    if(keep_order == 1) return(order.by.id...and.remove.it(merge(x=add.id.column.to.data(x),y=y,..., sort = FALSE)))
    if(keep_order == 2) return(order.by.id...and.remove.it(merge(x=x,y=add.id.column.to.data(y),..., sort = FALSE)))
    # if you didn't get "return" by now - issue a warning.
    warning("The function merge.with.order only accepts NULL/1/2 values for the keep_order variable")
  } else {return(merge(x=x,y=y,..., sort = sort))}
}

toPyr <- function(nctds) {
  nctds.pyr <- nctds
  nctds.pyr[nctds=='A'] <- 'T'
  nctds.pyr[nctds=='T'] <- 'A'
  nctds.pyr[nctds=='G'] <- 'C'
  nctds.pyr[nctds=='C'] <- 'G'
  as.character(nctds.pyr)
}
