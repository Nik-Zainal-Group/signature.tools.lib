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
