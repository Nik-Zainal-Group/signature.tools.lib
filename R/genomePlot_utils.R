processSubs <- function(subs) 
{
	# new style
	scatter.data  <-  calcIntermutDist(subs)
	scatter.data$mutType <- paste(scatter.data$ref_base_pyrimidine_context,'>',scatter.data$mutant_base_pyrimidine_context ,sep='')

	scatter.colors <- rep("gold", nrow(scatter.data)) #default color for multi-allelics
	scatter.colors[scatter.data$mutType=="C>A"] <- "royalblue"
	scatter.colors[scatter.data$mutType=="C>G"] <- "black"
	scatter.colors[scatter.data$mutType=="C>T"] <- "red"
	scatter.colors[scatter.data$mutType=="T>A"] <- "grey"
	scatter.colors[scatter.data$mutType=="T>C"] <- "green2"
	scatter.colors[scatter.data$mutType=="T>G"] <- "hotpink"

	result <- list()
	result$scatter.colors <- scatter.colors
	result$scatter.data <- scatter.data

	result

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


getMutTables <- function(myFile, onlyPASSED=FALSE, genome.v="hg19", genomeSeq, addContext=TRUE,mut.order) {
  # plots mutation-context for all variants in the vcf file
  # and separately for the variants that passed

if(genome.v=="hg19"){
  expected_chroms <- paste0(c(seq(1:22),"X","Y"))
}else if(genome.v=="hg38"){
  expected_chroms <- paste0("chr",c(seq(1:22),"X","Y"))
}else if(genome.v=="mm10"){
   expected_chroms <- paste0(c(seq(1:19),"X","Y")) 
}else if (genome.v=="canFam3"){
   expected_chroms <- paste0("chr",c(seq(1:38),"X"))  
}

# read only chr seqnames from VCF, not contigs
gr <- GenomicRanges::GRanges(GenomeInfoDb::seqinfo(genomeSeq))
if (genome.v=="hg19" || genome.v=="mm10") {
  GenomeInfoDb::seqlevels(gr) <- sub("chr", "", GenomeInfoDb::seqlevels(gr))
}
vcf_seqnames <- Rsamtools::headerTabix(myFile)$seqnames 

if(tools:::.BioC_version_associated_with_R_version()<3.5){
  gr <- GenomeInfoDb::keepSeqlevels(gr,intersect(vcf_seqnames,expected_chroms))
}else{
  gr <- GenomeInfoDb::keepSeqlevels(gr,intersect(vcf_seqnames,expected_chroms),pruning.mode = "coarse")
}

# load the vcf file
vcf_data <- VariantAnnotation::readVcf(myFile, genome.v, gr)

#browser()

#vcf_data <- keepSeqlevels(vcf_data, seqnames(genomeSeq))

#filters failed for each variant
rd <- SummarizedExperiment::rowData(vcf_data)

fs <- rd$FILTER
fs.passed <- (fs=='PASS')

if (onlyPASSED) {

    vcf_data <- vcf_data[fs.passed,]
    rd <- SummarizedExperiment::rowData(vcf_data)
    fs <- rd$FILTER
    fs.passed <- (fs=='PASS')

}

info.data <- VariantAnnotation::info(vcf_data)


rgs <- IRanges::ranges(vcf_data)
starts <- BiocGenerics::start(rgs)
ends <-  BiocGenerics::end(rgs)
chroms <- GenomeInfoDb::seqnames(vcf_data)
if (genome.v=="mm10"){
  chroms <- paste('chr',chroms,sep='')
}

wt <- as.character(rd$REF)
#mt <- as.character(unlist(rd$ALT))
mt <- IRanges::CharacterList(rd$ALT)
mt <- Biostrings::unstrsplit(mt, sep = ",")

barcode <- paste(chroms, '-',starts,'-', mt, sep='')

if (addContext) {
    bb <- as.character(BSgenome::getSeq(genomeSeq, chroms, start=starts-1, end=ends-1))
    ba <- as.character(BSgenome::getSeq(genomeSeq, chroms, start=starts+1, end=ends+1))
    wt.ref <- as.character(BSgenome::getSeq(genomeSeq, chroms, start=starts, end=ends))
    triplets <- as.character(BSgenome::getSeq(genomeSeq, chroms, start=starts-1, end=ends+1))

    # check the annotation
    if (sum(!wt.ref==wt)>0) {
        cat('wrong reference genome \n')
        browser()
    }

    mut.table <- data.frame(bbef=as.character(bb), wt=as.character(wt), mt=as.character(mt), baft=as.character(ba), stringsAsFactors=FALSE)


    mut.table$pyrwt <- as.character(mut.table$wt)
    mut.table$pyrmut <- as.character(mut.table$mt)
    mut.table$pyrbbef <- as.character(mut.table$bbef)
    mut.table$pyrbaft <- as.character(mut.table$baft)


    # the mutations originally not in pyramidine contex
    not.pyr <- ((wt=='G') | (wt=='A'))
    mut.table$pyrwt[not.pyr] <- as.character(toPyr(mut.table$wt[not.pyr]))
    mut.table$pyrmut[not.pyr] <- as.character(toPyr(mut.table$mt[not.pyr]))
    mut.table$pyrbbef[not.pyr] <- as.character(toPyr(mut.table$baft[not.pyr]))
    mut.table$pyrbaft[not.pyr] <- as.character(toPyr(mut.table$bbef[not.pyr]))

    all.hist <- generateHist(mut.table, normalise=FALSE,mut.order=mut.order)
    passed.hist <- generateHist(mut.table[fs.passed,], normalise=FALSE,mut.order=mut.order)
    # if (sum(info.data[,'SNP'])>0) {
    #     snps.hist <-  generateHist(mut.table[info.data[,'SNP'],], normalise=FALSE,mut.order=mut.order)
    #     names(snps.hist) <- mut.order
    # } else {
    #     snps.hist <- NULL
    # }
    # non.snps.hist <-  generateHist(mut.table[!info.data[,'SNP'],], normalise=FALSE,mut.order=mut.order)

    names(passed.hist) <- mut.order
    #names(non.snps.hist) <- mut.order

    muts <- data.frame(chroms=chroms, starts=starts, ends = ends, wt=wt, mt=mt, pyrwt=mut.table$pyrwt , pyrmut=mut.table$pyrmut, pass=fs.passed, barcode=barcode,
                   context=paste(mut.table$pyrbbef, '[',mut.table$pyrwt, '>',mut.table$pyrmut , ']', mut.table$pyrbaft,sep=''),
                   tumor.freq=rep(NA,length(chroms)), normal.freq=rep(NA,length(chroms)),
                   tumor.reads=rep(NA,length(chroms)), normal.reads=rep(NA,length(chroms)),
                   tumor.depth=rep(NA,length(chroms)), normal.depth=rep(NA,length(chroms)),
                                                           isSnp =rep(NA,length(chroms))
                   )
} else {

    mut.table <- NULL

    muts <- data.frame(chroms=chroms, starts=starts, ends = ends, wt=wt, mt=mt,
                       barcode=barcode,
                       tumor.freq=rep(NA,length(chroms)), normal.freq=rep(NA,length(chroms)),
                       tumor.reads=rep(NA,length(chroms)), normal.reads=rep(NA,length(chroms)),
                       tumor.depth=rep(NA,length(chroms)), normal.depth=rep(NA,length(chroms)),
                       isSnp =rep(NA,length(chroms)), filters=fs
                       )

    all.hist <- NULL
    passed.hist <- NULL
    mut.table <- NULL
    #snps.hist <- NULL
    #non.snps.hist <- NULL

}

# calculate tumour frequency

# calculate normal frequency
geno.data <- VariantAnnotation::geno(vcf_data)

if ('FAZ' %in% names(geno.data)) {
    muts$tumor.depth <- geno.data[['FAZ']][,'TUMOUR'] + geno.data[['RAZ']][,'TUMOUR'] +
        geno.data[['FCZ']][,'TUMOUR'] + geno.data[['RCZ']][,'TUMOUR'] +
            geno.data[['FGZ']][,'TUMOUR'] + geno.data[['RGZ']][,'TUMOUR'] +
                geno.data[['FTZ']][,'TUMOUR'] + geno.data[['RTZ']][,'TUMOUR']
    muts$normal.depth <- geno.data[['FAZ']][,'NORMAL'] + geno.data[['RAZ']][,'NORMAL'] +
        geno.data[['FCZ']][,'NORMAL'] + geno.data[['RCZ']][,'NORMAL'] +
            geno.data[['FGZ']][,'NORMAL'] + geno.data[['RGZ']][,'NORMAL'] +
                geno.data[['FTZ']][,'NORMAL'] + geno.data[['RTZ']][,'NORMAL']

    is.A <- (mt=='A')
    muts$tumor.reads[is.A] <- ((geno.data[['FAZ']][,'TUMOUR'][is.A] +  geno.data[['RAZ']][,'TUMOUR'][is.A]))
    muts$normal.reads[is.A] <- ((geno.data[['FAZ']][,'NORMAL'][is.A] +  geno.data[['RAZ']][,'NORMAL'][is.A]))
    muts$tumor.freq[is.A] <- ((geno.data[['FAZ']][,'TUMOUR'][is.A] +  geno.data[['RAZ']][,'TUMOUR'][is.A]))/ muts$tumor.depth[is.A]
    muts$normal.freq[is.A] <- ((geno.data[['FAZ']][,'NORMAL'][is.A] +  geno.data[['RAZ']][,'NORMAL'][is.A]))/ muts$normal.depth[is.A]

    is.C <- (mt=='C')
    muts$tumor.reads[is.C] <- ((geno.data[['FCZ']][,'TUMOUR'][is.C] +  geno.data[['RCZ']][,'TUMOUR'][is.C]))
    muts$normal.reads[is.C] <- ((geno.data[['FCZ']][,'NORMAL'][is.C] +  geno.data[['RCZ']][,'NORMAL'][is.C]))
    muts$tumor.freq[is.C] <- ((geno.data[['FCZ']][,'TUMOUR'][is.C] +  geno.data[['RCZ']][,'TUMOUR'][is.C]))/ muts$tumor.depth[is.C]
    muts$normal.freq[is.C] <- ((geno.data[['FCZ']][,'NORMAL'][is.C] +  geno.data[['RCZ']][,'NORMAL'][is.C]))/ muts$normal.depth[is.C]

    is.G <- (mt=='G')
    muts$tumor.reads[is.G] <- ((geno.data[['FGZ']][,'TUMOUR'][is.G] +  geno.data[['RGZ']][,'TUMOUR'][is.G]))
    muts$normal.reads[is.G] <- ((geno.data[['FGZ']][,'NORMAL'][is.G] +  geno.data[['RGZ']][,'NORMAL'][is.G]))
    muts$tumor.freq[is.G] <- ((geno.data[['FGZ']][,'TUMOUR'][is.G] +  geno.data[['RGZ']][,'TUMOUR'][is.G]))/ muts$tumor.depth[is.G]
    muts$normal.freq[is.G] <- ((geno.data[['FGZ']][,'NORMAL'][is.G] +  geno.data[['RGZ']][,'NORMAL'][is.G]))/ muts$normal.depth[is.G]

    is.T <- (mt=='T')
    muts$tumor.reads[is.T] <- ((geno.data[['FTZ']][,'TUMOUR'][is.T] +  geno.data[['RTZ']][,'TUMOUR'][is.T]))
    muts$normal.reads[is.T] <- ((geno.data[['FTZ']][,'NORMAL'][is.T] +  geno.data[['RTZ']][,'NORMAL'][is.T]))
    muts$tumor.freq[is.T] <- ((geno.data[['FTZ']][,'TUMOUR'][is.T] +  geno.data[['RTZ']][,'TUMOUR'][is.T]))/ muts$tumor.depth[is.T]
    muts$normal.freq[is.T] <- ((geno.data[['FTZ']][,'NORMAL'][is.T] +  geno.data[['RTZ']][,'NORMAL'][is.T]))/ muts$normal.depth[is.T]
}

# whether the sub was found to be a SNP
#muts$isSnp <- info.data[,'SNP']

result<- list()
result$all.hist <- all.hist
result$passed.hist <- passed.hist
result$muts <- muts
result$mut.table <- mut.table
#result$snps.hist <- snps.hist
#result$non.snps.hist <- non.snps.hist
result

}


getMutTablesTab <- function(SUBS.PATH, onlyPASSED=FALSE, genomeSeq=Hsapiens, addContext=TRUE,mut.order) {
	# plots mutation-context for all variants in the vcf file
	# and separately for the variants that passed

    # load the tab file
    #Ensure the wt and mt columns (V7 & V8) are read in as characters.
    subs <- read.table(SUBS.PATH, header=FALSE, sep='\t', quote = "", stringsAsFactors=FALSE, colClasses=c(V7="character",V8="character"))
    subs$chr <- as.character(subs$V5)
    subs$chr.chr<- paste('chr', subs$chr, sep='')
    subs$position <- subs$V6
    subs$wt <- subs$V7
    subs$mt <- subs$V8
    subs$ref_base_pyrimidine_context <- subs$wt
    subs$mutant_base_pyrimidine_context <- subs$mt
    noPyrBases <- (subs$wt=='G') | (subs$wt=='A')
    subs$ref_base_pyrimidine_context[noPyrBases ] <- toPyr(subs$wt[noPyrBases])
    subs$mutant_base_pyrimidine_context[noPyrBases ]  <- toPyr(subs$mt[noPyrBases])

    fs.passed<-   subs$V10=='PASS'
    if (onlyPASSED) {
        subs <- subset(subs, V10=='PASS')
    }

                                        # the mutations originally not in pyramidine contex
    not.pyr <- ((subs$wt=='G') | (subs$wt=='A'))
    subs$pyrwt <- as.character(subs$wt)
    subs$pyrmut <- as.character(subs$mt)
    subs$pyrwt[not.pyr] <- as.character(toPyr(subs$wt[not.pyr]))
    subs$pyrmut[not.pyr] <- as.character(toPyr(subs$mt[not.pyr]))

    barcode <- paste(subs$chr, '-',subs$position ,'-', subs$mt, sep='')

    muts <- data.frame(chroms=subs$chr.chr, starts=subs$position, ends = subs$position, wt=subs$wt, mt=subs$mt, pyrwt=subs$pyrwt , pyrmut=subs$pyrmut, pass=subs$V10, barcode=barcode
                   )


    result<- list()
    result$muts <- muts
######

    if (addContext) {

        subs$bb <- as.character(getSeq(genomeSeq, subs$chr.chr, start=subs$position-1, end=subs$position-1))
        subs$ba <- as.character(getSeq(genomeSeq, subs$chr.chr, start=subs$position+1, end=subs$position+1))
        subs$wt.ref <- as.character(getSeq(genomeSeq, subs$chr.chr, start=subs$position, end=subs$position))
        subs$triplets <- as.character(getSeq(genomeSeq, subs$chr.chr, start=subs$position-1, end=subs$position+1))

                                        # table of mutations
        mut.table <- data.frame(bbef=as.character(subs$bb ), wt=as.character(subs$wt), mt=as.character(subs$mt), baft=as.character(subs$ba ))
        mut.table$pyrwt <- mut.table$wt
        mut.table$pyrmut <- mut.table$mt
        mut.table$pyrbbef <- mut.table$bbef
        mut.table$pyrbaft <- mut.table$baft

                                        # the mutations originally not in pyramidine contex
        mut.table$pyrwt[not.pyr] <- as.character(toPyr(mut.table$wt[not.pyr]))
        mut.table$pyrmut[not.pyr] <- as.character(toPyr(mut.table$mt[not.pyr]))
        mut.table$pyrbbef[not.pyr] <- as.character(toPyr(mut.table$baft[not.pyr]))
        mut.table$pyrbaft[not.pyr] <- as.character(toPyr(mut.table$bbef[not.pyr]))

        all.hist <- generateHist(mut.table, normalise=FALSE,mut.order=mut.order)
        passed.hist <- generateHist(mut.table[fs.passed,], normalise=FALSE,mut.order=mut.order)

        names(passed.hist) <- mut.order

        muts$context <- paste(mut.table$pyrbbef, '[',mut.table$pyrwt, '>',mut.table$pyrmut , ']', mut.table$pyrbaft,sep='')

        result$all.hist <- all.hist
        result$passed.hist <- passed.hist
        result$muts <- muts
    }

    result

}


# function for importing the result of ASCAT or ASCAT NGS.
# arguments
# FILE.CN - full path the the ascat output (extension .ascat.summary.csv or ascat_ngs.summary.csv)
read.ascat <- function(FILE.CN) 
{
    #first row
    firstline <- tryCatch(read.table(text=gsub(",", "\t",readLines(FILE.CN,n = 1)),header=FALSE, sep='\t',nrows = 1,stringsAsFactors = FALSE), error=function(e) NULL)
    if(typeof(firstline[1,1])=="integer"){
      #no header used
  	  cv.data <- tryCatch(read.table(text = gsub(",", "\t", readLines(FILE.CN)),header=FALSE, sep='\t',stringsAsFactors = FALSE), error=function(e) NULL) # ASCAT
  	  names(cv.data) <- c('seg_no', 'Chromosome', 'chromStart', 'chromEnd', 'total.copy.number.inNormal', 'minor.copy.number.inNormal', 'total.copy.number.inTumour', 'minor.copy.number.inTumour')
    }else{
      #assume there is an appropriate header
      cv.data <- tryCatch(read.table(text = gsub(",", "\t", readLines(FILE.CN)),header=TRUE, sep='\t',stringsAsFactors = FALSE), error=function(e) NULL) # CN
    }
    cat( paste( dim(cv.data)[1], ' copy-number segments \n'))

    if (!is.null(cv.data) && !inherits(cv.data, "try-error")) {
        cv.data$seg_no <- NULL
        cv.data$Chromosome <- as.character(cv.data$Chromosome)
        #cv.data$Chromosome[cv.data$Chromosome=='23'] <- 'X'
        #cv.data$Chromosome[cv.data$Chromosome=='24'] <- 'Y'    
        cv.data$major.copy.number.inTumour <- cv.data$total.copy.number.inTumour - cv.data$minor.copy.number.inTumour
        cv.data$major.copy.number.inTumour.temp <- pmax(cv.data$major.copy.number.inTumour, cv.data$minor.copy.number.inTumour)
        cv.data$minor.copy.number.inTumour.temp <- pmin(cv.data$major.copy.number.inTumour, cv.data$minor.copy.number.inTumour)
        
        cv.data$major.copy.number.inTumour <- cv.data$major.copy.number.inTumour.temp
        cv.data$minor.copy.number.inTumour <- cv.data$minor.copy.number.inTumour.temp
        return(cv.data)
    } else {
        return(data.frame())
    }    

}


read.brassII <- function(FILE.REARR) 
{
    rearrs <- tryCatch(suppressWarnings(read.table(FILE.REARR, header=FALSE, sep='\t')), error=function(e) NULL)
    if (!inherits(rearrs, "try-error")) {
        rearrs$V1 <- as.character(rearrs$V1 )
        rearrs$V5 <- as.character(rearrs$V5 )
        #rearrs$V1 [rearrs$V1=='23'] <- 'X'
        #rearrs$V5[rearrs$V5=='24'] <- 'Y'
        cat(paste( dim(rearrs)[1], ' rearrs \n'))
        
        pf <- rep(4,nrow(rearrs))
        f.32 <- rearrs$V1!=rearrs$V5; pf[f.32] <- 32 # translocations
        f.1 <- rearrs$V1==rearrs$V5 & rearrs$V2=='+' & rearrs$V6=='-'; pf[f.1] <- 1 # inversion +-
        f.8 <- rearrs$V1==rearrs$V5 & rearrs$V2=='-' & rearrs$V6=='+'; pf[f.8] <- 8 # inversion -+
        f.2 <- rearrs$V1==rearrs$V5 & rearrs$V2=='+' & rearrs$V6=='+'; pf[f.2] <- 2 # deletion
        rearrs$pf <- pf
        
        rearrs.formatted <- data.frame(Chromosome=rearrs$V1, chromStart=rearrs$V3, chromEnd=rearrs$V3, Chromosome.1=rearrs$V5, chromStart.1=rearrs$V7, chromEnd.1=rearrs$V7, pf=rearrs$pf)
        return(rearrs.formatted)
    } else {
        return(data.frame())
    }   
}


read.brass.bedpe <- function(FILE.REARR, onlyAssembled = TRUE) 
{
    rearrs <- tryCatch(read.table(gzfile(FILE.REARR), header=TRUE, sep='\t',  comment.char = ''), error=function(e) NULL)

    if (!is.null(rearrs) && !inherits(rearrs, "try-error")) {
        names(rearrs)[1] <- 'chr1'
        names(rearrs)[2] <- 'start1'
        names(rearrs)[3] <- 'end1'
        names(rearrs)[4] <- 'chr2'
        names(rearrs)[5] <- 'start2'
        names(rearrs)[6] <- 'end2'

        rearrs$chr1 <- as.character(rearrs$chr1)
        rearrs$chr2 <- as.character(rearrs$chr2)


        pf <- rep(4,nrow(rearrs))
        
        if ("TYPE" %in% names(rearrs)) 
        {
            f.32 <- rearrs$TYPE=="BND"; pf[f.32] <- 32 # translocations
            f.1  <- rearrs$TYPE=="INV"; pf[f.1] <- 1 # inversion +-
            f.2  <- rearrs$TYPE=="DEL"; pf[f.2] <- 2 # deletion ++
        } 
        else if ("svclass" %in% names(rearrs)) 
        {
          f.32 <- rearrs$svclass=="translocation"; pf[f.32] <- 32 # translocations
          f.1  <- rearrs$svclass=="inversion"; pf[f.1] <- 1 # inversion +-
          f.2  <- rearrs$svclass=="deletion"; pf[f.2] <- 2 # deletion ++
        } 
        else 
        {
            f.32 <- rearrs$chr1!=rearrs$chr2; pf[f.32] <- 32 # translocations
            f.1  <- rearrs$chr1==rearrs$chr2 & rearrs$strand1=='+' & rearrs$strand2=='-'; pf[f.1] <- 1 # inversion +-
            f.8  <- rearrs$chr1==rearrs$chr2 & rearrs$strand1=='-' & rearrs$strand2=='+'; pf[f.8] <- 8 # inversion -+
            f.2  <- rearrs$chr1==rearrs$chr2 & rearrs$strand1=='+' & rearrs$strand2=='+'; pf[f.2] <- 2 # deletion ++
        }

        rearrs$pf <- pf

        rearrs.formatted <- data.frame(Chromosome=rearrs$chr1,
                                       chromStart=apply(rearrs[,c('start1','end1')],1, min),
                                       chromEnd=apply(rearrs[,c('start1','end1')],1, max),
                                       Chromosome.1=rearrs$chr2,
                                       chromStart.1=apply(rearrs[,c('start2','end2')],1, min),
                                       chromEnd.1=apply(rearrs[,c('start2','end2')],1, max),
                                       pf=rearrs$pf)

        if (onlyAssembled) {
            if (nrow(rearrs.formatted)>0) {
                rearrs.formatted  <- rearrs.formatted[rearrs$assembly_score!='_', ] # only include the rearrangements that have an assembly score
            }
        }
        cat(paste( dim(rearrs.formatted)[1], ' rearrs \n'))
        return(rearrs.formatted)
    } else {
        return(data.frame())
    }    
}


read.brassII <- function(FILE.REARR) {


    rearrs <- tryCatch(suppressWarnings(read.table(FILE.REARR, header=FALSE, sep='\t')), error=function(e) NULL)
    if (!inherits(rearrs, "try-error")) {
        rearrs$V1 <- as.character(rearrs$V1 )
        rearrs$V5 <- as.character(rearrs$V5 )
        #rearrs$V1 [rearrs$V1=='23'] <- 'X'
        #rearrs$V5[rearrs$V5=='24'] <- 'Y'
        cat(paste( dim(rearrs)[1], ' rearrs \n'))
        
        pf <- rep(4,nrow(rearrs))
        f.32 <- rearrs$V1!=rearrs$V5; pf[f.32] <- 32 # translocations
        f.1 <- rearrs$V1==rearrs$V5 & rearrs$V2=='+' & rearrs$V6=='-'; pf[f.1] <- 1 # inversion +-
        f.8 <- rearrs$V1==rearrs$V5 & rearrs$V2=='-' & rearrs$V6=='+'; pf[f.8] <- 8 # inversion -+
        f.2 <- rearrs$V1==rearrs$V5 & rearrs$V2=='+' & rearrs$V6=='+'; pf[f.2] <- 2 # deletion
        rearrs$pf <- pf
        
        rearrs.formatted <- data.frame(Chromosome=rearrs$V1, chromStart=rearrs$V3, chromEnd=rearrs$V3, Chromosome.1=rearrs$V5, chromStart.1=rearrs$V7, chromEnd.1=rearrs$V7, pf=rearrs$pf)
        return(rearrs.formatted)
    } else {
        return(data.frame())
    }
    
    
}


