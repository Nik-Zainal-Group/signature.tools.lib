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
