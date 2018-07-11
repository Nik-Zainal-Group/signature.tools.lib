### Nicolai Juul Birkbak, njuul@cbs.dtu.dk
### Functions to Serena/Johan, 2014-07-02
## To calculate NtAI, LST, HRD
# These functions generally rely on ASCAT or similarly processed copy number data, in a matrix format with at least the follwing columns:
# SampleID	Chromosome	Start	End	nProbes	totalCN	nA	nB	Ploidy	AberrantCellFraction ## more columns may be added. nA and nB stands for allelic copy number (ie., CN of allele A, and CN of allele B)

########## VERY important!!!!! ##########
# Note that the current implementation of calc.lst requires a different format for the chromosome info file!! I plan to update the other functions to take the same format, but haven't gotten to it yet.
# Hence, use the "chrominfo" file for all functions except calc.lst. For calc.lst, use chrominfo.snp6. 
#########################################


# Feel free to use as you please. If you find any bugs, please get back to me!


################
#Number of telomeric AI, updated and improved to also account for ploidy - 2013-03-11
################
#' Number of Telomeric Allelic Imbalances
#' 
#' This function computes the number of Telomeric AI (NtAI). Implementation of Nicolai Juul Birkbak, njuul@cbs.dtu.dk (2014-07-02).
#' This functions generally rely on ASCAT or similarly processed copy number data, in a matrix format with at least the follwing columns in this exact order:
#' "SampleID", "Chromosome", "Start", "End", "nProbes", "totalCN", "nA", "nB", "Ploidy" and	"AberrantCellFraction".
#' NOTE: currently the chrominfo data refers to hg19, while hg38 is not yet supported.
#' 
#' @param seg input data frame. Segmented output in the form of an ASCAT out matrix, which also includes ploidy and contamination. Total CN in column 6, nA in column 7, nB in column 8, ploidy in column 9, and contamination in column 10. Requires the following columns in the given order: "SampleID", "Chromosome", "Start", "End", "nProbes", "totalCN", "nA", "nB", "Ploidy" and	"AberrantCellFraction". 
#' @param chrominfo data frame that contains necessary chromosome information. No need to specify this for human chromosomes, it will be loaded internally. A 3 column matrix with information about the chromosomes: chromosome name, chromosome length, centromere location.
#' @param min.size minimum size of segments
#' @param min.probes minimum number of probes in segments (I use 500 for the 900,000 probe SNP6, then scale down)
#' @param cont contamination threshold. By default set at 0 to ignore contamination
#' @param check.names check and potentially fix if there are duplicated samples. Any duplicates shall be re-named. Works by combining sample names with ploidy and SNP file name (which for this to work must be in column 11). Then looks for individual sample names which matches multiple original SNP files.
#' @param shrink joins segments of identical allelic copy number
#' @export
#' @author Nicolai Juul Birkbak, njuul@cbs.dtu.dk
#' @references Birkbak, N. J., Wang, Z. C., Kim, J.-Y., Eklund, A. C., Li, Q., Tian, R., ... Richardson, A. L. (2012). Telomeric Allelic Imbalance Indicates Defective DNA Repair and Sensitivity to DNA-Damaging Agents. Cancer Discovery, 2(4), 366–375. https://doi.org/10.1158/2159-8290.CD-11-0206
#' @examples 
#' res <- calc.ai(seg)
calc.ai <- function(seg, chrominfo=chrominfo, min.size=0, min.probes=500, cont = 0, check.names=FALSE, ploidyByChromosome=TRUE, return.loc=FALSE, shrink=TRUE){
	#edit 2014-01-17, added the ability to return location of NtAI's (return.loc)
	#edit 2014-02-12, added the "shrink" option, which joins segments with same allelic copy number. These may be brought together if segments in between are filtered due to some settings, like minimum probe numbers, or size
	#edit 2014-02-25, check.names is now a separate function
	#edit 2014-05-27, disallow ploidy of zero
	
	# seg = segmented output in the form of an ASCAT out matrix, which also includes ploidy and contamination. Total CN in column 6, nA in column 7, nB in column 8, ploidy in column 9, and contamination in column 10
	# chrominfo = a 3 column matrix with information about the chromosomes: chromosome name, chromosome length, centromere location.
	# min.size = minimum size of segments
	# min.probes = minimum number of probes in segments (I use 500 for the 900,000 probe SNP6, then scale down)
	# type = should the algorithm test for AI or LOH?
	# cont = contamination threshold. By default set at 0 to ignore contamination
	# check.names = check and potentially fix if there are duplicated samples. Any duplicates shall be re-named. Works by combining sample names with ploidy and SNP file name (which for this to work must be in column 11). Then looks for individual sample names which matches multiple original SNP files.
	# shrink = joins segments of identical allelic copy number
	
	if(ploidyByChromosome){cat("Determining chromosome-specific ploidy by major copy number fraction\n")}
	if(!ploidyByChromosome){cat("Determining sample ploidy by major copy number fraction over-all\n")}
	seg[,1] <- as.character(seg[,1])
	if(check.names){
		seg <- check.names.fn(seg)
	}
	# check and potentially fix if nB is always the smaller column
	if(!all(seg[,8] <= seg[,7]) ){
		cat("Warning!! nB  not always <= nA!!  -- Correcting for internal use (only!)\n") # In case ASCAT people change the algorithm
		tmp <- seg
		seg[tmp[,8] > tmp[,7],7]  <- tmp[tmp[,8] > tmp[,7],8]
		seg[tmp[,8] > tmp[,7],8]  <- tmp[tmp[,8] > tmp[,7],7]
	}
	# remove segments smaller min.size and min.probes, and with too much contamination ##### NOTE: for backward-compatibility, this might need to be moved to after calling of telomeric segments
	seg <- seg[seg[,5] >= min.probes,]
	seg <- seg[seg[,4]- seg[,3] >= min.size,]
	seg <- seg[seg[,10] >= cont,]
	samples <- as.character(unique(seg[,1]))
	if(shrink){
		new.seg <- seg[1,]
		for(j in samples){
			sample.seg <- seg[seg[,1] %in% j,]
			new.sample.seg <- seg[1,]
			for(i in unique(sample.seg[,2])){
				sample.chrom.seg <- sample.seg[sample.seg[,2] %in% i,,drop=F]
				sample.chrom.seg <- shrink.seg.ai(sample.chrom.seg)
				new.sample.seg <- rbind(new.sample.seg, sample.chrom.seg)
			}
			new.seg <- rbind(new.seg, new.sample.seg[-1,])		
		}
		seg <- new.seg[-1,]
	}			
	#Add a column to call AI
	AI <- rep(NA, nrow(seg))
	seg <- cbind(seg, AI)
	samples <- as.character(unique(seg[,1]))
	ascat.ploidy <- setNames(seg[!duplicated(seg[,1]),9], seg[!duplicated(seg[,1]),1])
	for(j in samples){
		sample.seg <- seg[seg[,1] %in% j,]
		if(!ploidyByChromosome){
			ploidy <- vector()
			for(k in unique(sample.seg[,6])){
				tmp <- sample.seg[sample.seg[,6] %in% k,]
				ploidy <- c(ploidy, setNames(sum(tmp[,4]-tmp[,3]), k))
			}
			ploidy <- as.numeric(names(ploidy[order(ploidy,decreasing=T)]))[1]		
			sample.seg[,9] <- ploidy # update "ploidy" column, so the new calculated value can be returned
			# add a columnm to define AI, with codes for telomeric/interstitial/whole chromosome. 1= telomeric, 2= interstitial, 3 = whole chromosome
			if(ploidy %in% c(1,seq(2, 200,by=2))){
				sample.seg[,'AI'] <- c(0,2)[match(sample.seg[,7] == sample.seg[,8], c('TRUE', 'FALSE'))]
			}		
			if(!ploidy %in%  c(1,seq(2, 200,by=2))){
				sample.seg[,'AI'] <- c(0,2)[match(sample.seg[,7] + sample.seg[,8] == ploidy & sample.seg[,7] != ploidy, c('TRUE', 'FALSE'))]
			}
		}
		new.sample.seg<- sample.seg[1,]
		for(i in unique(sample.seg[,2])){
			sample.chrom.seg <- sample.seg[sample.seg[,2] %in% i,,drop=F]
			if(nrow(sample.chrom.seg) == 0){ next}
			if(ploidyByChromosome){
				ploidy <- vector()
				for(k in unique(sample.seg[,6])){
					tmp <- sample.chrom.seg[sample.chrom.seg[,6] %in% k,]
					ploidy <- c(ploidy, setNames(sum(tmp[,4]-tmp[,3]), k))
					#Remove any ploidy calls of zero
					ploidy <- ploidy[!names(ploidy) %in% 0]
				}
				ploidy <- as.numeric(names(ploidy[order(ploidy,decreasing=T)]))[1]		
				sample.chrom.seg[,9] <- ploidy # update "ploidy" column, so the new calculated value can be returned
				if(ploidy %in% c(1,seq(2, 200,by=2))){
					sample.chrom.seg[,'AI'] <- c(0,2)[match(sample.chrom.seg[,7] == sample.chrom.seg[,8], c('TRUE', 'FALSE'))]
				}		
				if(!ploidy %in%  c(1,seq(2, 200,by=2))){
					sample.chrom.seg[,'AI'] <- c(0,2)[match(sample.chrom.seg[,7] + sample.chrom.seg[,8] == ploidy & sample.chrom.seg[,8] != 0, c('TRUE', 'FALSE'))]
				}
				sample.seg[sample.seg[,2] %in% i,9] <-ploidy
				sample.seg[sample.seg[,2] %in% i,'AI'] <-sample.chrom.seg[,'AI']
			}
			if(class(chrominfo) == 'logical'){			# By logical, we assume that chrominfo == FALSE, hence we here proceed without checking for the centromere (useful for non-human samples)
				if(sample.chrom.seg[1,'AI'] == 2 & nrow(sample.chrom.seg) != 1){
					sample.seg[sample.seg[,2]==i,'AI'][1] <- 1
				}
				if(sample.chrom.seg[nrow(sample.chrom.seg),'AI'] == 2 & nrow(sample.chrom.seg) != 1){
					sample.seg[sample.seg[,2]==i,'AI'][nrow(sample.seg[sample.seg[,2]==i,])] <- 1
				}
			}
			if(class(chrominfo) != 'logical'){		# Here we consider the centromere		
				if(sample.chrom.seg[1,'AI'] == 2 & nrow(sample.chrom.seg) != 1 & sample.chrom.seg[1,4] < (chrominfo[i,3]*1000)){
					sample.seg[sample.seg[,2]==i,'AI'][1] <- 1
				}
				if(sample.chrom.seg[nrow(sample.chrom.seg),'AI'] == 2 & nrow(sample.chrom.seg) != 1 & sample.chrom.seg[nrow(sample.chrom.seg),3] > (chrominfo[i,3]*1000)){
					sample.seg[sample.seg[,2]==i,'AI'][nrow(sample.seg[sample.seg[,2]==i,])] <- 1
				}
			}
			if(nrow(sample.seg[sample.seg[,2]==i,]) == 1 & sample.seg[sample.seg[,2]==i,'AI'][1] != 0){
				sample.seg[sample.seg[,2]==i,'AI'][1] <- 3
			}
		}
		seg[seg[,1] %in% j,] <- sample.seg
	}	
	samples <- as.character(unique(seg[,1]))
	#0 = no AI, 1=telomeric AI, 2=interstitial AI, 3= whole chromosome AI
	no.events <- matrix(0, nrow=length(samples), ncol=12)
	rownames(no.events) <- samples
	colnames(no.events) <- c("Telomeric AI", "Mean size", "Interstitial AI", "Mean Size", "Whole chr AI", "Telomeric LOH",  "Mean size", "Interstitial LOH", "Mean Size", "Whole chr LOH", "Ploidy", "Aberrant cell fraction")
	for(j in samples){
		sample.seg <- seg[seg[,1] %in% j,]
		no.events[j,1] <- nrow(sample.seg[sample.seg[,'AI'] == 1,])
		no.events[j,2] <- mean(sample.seg[sample.seg[,'AI'] == 1,4] - sample.seg[sample.seg[,'AI'] == 1,3])
		no.events[j,3] <- nrow(sample.seg[sample.seg[,'AI'] == 2,])
		no.events[j,4] <- mean(sample.seg[sample.seg[,'AI'] == 2,4] - sample.seg[sample.seg[,'AI'] == 2,3])
		no.events[j,5] <- nrow(sample.seg[sample.seg[,'AI'] == 3,])
		no.events[j,11] <- ascat.ploidy[j]
		no.events[j,12] <- unique(sample.seg[,10]) # aberrant cell fraction
		#Here we restrict ourselves to real LOH
		sample.seg <- sample.seg[sample.seg[,8] == 0,]
		no.events[j,6] <- nrow(sample.seg[sample.seg[,'AI'] == 1,])
		no.events[j,7] <- mean(sample.seg[sample.seg[,'AI'] == 1,4] - sample.seg[sample.seg[,'AI'] == 1,3])
		no.events[j,8] <- nrow(sample.seg[sample.seg[,'AI'] == 2,])
		no.events[j,9] <- mean(sample.seg[sample.seg[,'AI'] == 2,4] - sample.seg[sample.seg[,'AI'] == 2,3])
		no.events[j,10] <- nrow(sample.seg[sample.seg[,'AI'] == 3,])
	}
	if(return.loc){
		return(seg)
	} else {
		return(no.events)
	}
}





##################
#This function implements Popova's LST measure as described (Popova 2012, Cancer research)
#Popova's cutoffs: 15 LSTs in near-diploid, 20 in near-tetraploid
##################
#' Large-scale state transitions
#' 
#' This function implements Popova's LST measure as described (Popova 2012, Cancer research)
#' Popova's cutoffs: 15 LSTs in near-diploid, 20 in near-tetraploid.
#' This functions generally rely on ASCAT or similarly processed copy number data, in a matrix format with at least the follwing columns in this exact order:
#' "SampleID", "Chromosome", "Start", "End", "nProbes", "totalCN", "nA", "nB", "Ploidy" and	"AberrantCellFraction".
#' NOTE: currently the chrominfo data refers to hg19, while hg38 is not yet supported.
#' 
#' @param seg seg must be an ASCAT output object, in DNAcopy format. Requires the following columns in the given order: "SampleID", "Chromosome", "Start", "End", "nProbes", "totalCN", "nA", "nB", "Ploidy" and	"AberrantCellFraction". Column "nProbes" is not used here and can be set to NA.
#' @param nA is the column where copy number of A allele is found. This needs to be set to 7, because of the columns are referenced by position according to the required column order of seg.
#' @export
#' @author Nicolai Juul Birkbak, njuul@cbs.dtu.dk
#' @references Popova, T., Manié, E., Rieunier, G., Caux-Moncoutier, V., Tirapo, C., Dubois, T., ... Stern, M. H. (2012). Ploidy and large-scale genomic instability consistently identify basal-like breast carcinomas with BRCA1/2 inactivation. Cancer Research, 72(21), 5454–5462. https://doi.org/10.1158/0008-5472.CAN-12-1470
#' @examples
#' res <- calc.lst(seg)
calc.lst <- function(seg, chrominfo=chrominfo.snp6,nA=7, check.names=FALSE, return.loc=FALSE, chr.arm='no'){
	#Seg must be an ASCAT output object, in DNAcopy format. 
	#nA is the column where copy number of A allele is found
	#Edit 20140216: added return.loc option, it will return the location and arm of the LST sites
	#Edit 20140325: added centromere locations based on SNP6 array: load('~/Desktop/DFProjects/GenomeData/chrominfo.snp6.RData')
	#Edit 20140528: Added the option to use chromosome arms defined during segmentation. The option must give a column that holds the chromosome arm information. 
	
	seg[,1] <- as.character(seg[,1])
	if(check.names){
		seg <- check.names.fn(seg)
	}
	if(! all(seg[,8] <= seg[,7]) ){
		cat("Warning!! nB  not always <= nA!!  -- Correcting for internal use (only!)\n") # In case ASCAT people change the algorithm  ### They did!!
		tmp <- seg
		seg[tmp[,8] > tmp[,7],7]  <- tmp[tmp[,8] > tmp[,7],8]
		seg[tmp[,8] > tmp[,7],8]  <- tmp[tmp[,8] > tmp[,7],7]
	}
	nB <- nA+1
	samples <- unique(seg[,1])
	output <- setNames(rep(0,length(samples)), samples)
	if(return.loc) {
			out.seg <- matrix(0,0,10)
			colnames(out.seg) <- c(colnames(seg)[1:8],'LST breakpoint', 'chr. arm')
	}
	for(j in samples){
		sample.seg <- seg[seg[,1] %in% j,]
		sample.lst <- c()
		for(i in unique(sample.seg[,2])){
			sample.chrom.seg <- sample.seg[sample.seg[,2] %in% i,,drop=F]
			if(chr.arm !='no'){
				p.max <- if(any(sample.chrom.seg[,chr.arm] == 'p')){max(sample.chrom.seg[sample.chrom.seg[,chr.arm] == 'p',4])}
				q.min <- min(sample.chrom.seg[sample.chrom.seg[,chr.arm] == 'q',3])
			}
			sample.chrom.seg <- sample.chrom.seg[(sample.chrom.seg[,4] - sample.chrom.seg[,3]) > 3e6,] # remove breaks shorter than 3 mb 
			if(nrow(sample.chrom.seg) < 2) {next}
			sample.chrom.seg.new <- sample.chrom.seg
#			sample.chrom.seg.new <- sample.chrom.seg[1,,drop=F]
#			count <- 1
#			for(k in 2:nrow(sample.chrom.seg)){
#				if(!(sample.chrom.seg[k,7] == sample.chrom.seg.new[count,7] & sample.chrom.seg[k,8] == sample.chrom.seg.new[count,8])){
#					sample.chrom.seg.new <- rbind(sample.chrom.seg.new, sample.chrom.seg[k,,drop=F])
#					sample.chrom.seg.new[count,4] <- sample.chrom.seg.new[count+1,3]
#					count<- count+1
#				}
#			}
#			if(sample.chrom.seg[nrow(sample.chrom.seg),4] != sample.chrom.seg.new[nrow(sample.chrom.seg.new),4]){
#				sample.chrom.seg.new[nrow(sample.chrom.seg.new),4] <- sample.chrom.seg[nrow(sample.chrom.seg),4]
#			}
#			sample.chrom.seg <- cbind(sample.chrom.seg.new[,1:8], c(0,1)[match((sample.chrom.seg.new[,4]-sample.chrom.seg.new[,3]) >= 10e6, c('FALSE','TRUE'))])
			if(chr.arm == 'no'){
				p.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,3] <= chrominfo[i,4],] # split into chromosome arms
				q.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,4] >= chrominfo[i,5],]
				q.arm<- shrink.seg.ai(q.arm)
				p.arm<- shrink.seg.ai(p.arm)
				p.arm[nrow(p.arm),4] <- chrominfo[i,4]
				q.arm[1,3] <- chrominfo[i,5]
			}
			if(chr.arm != 'no'){
				q.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,chr.arm] == 'q',,drop=F]
				q.arm<- shrink.seg.ai(q.arm)
				q.arm[1,3] <- q.min
				if(any(sample.chrom.seg.new[,chr.arm] == 'p')){
					p.arm <- sample.chrom.seg.new[sample.chrom.seg.new[,chr.arm] == 'p',,drop=F] # split into chromosome arms
					p.arm<- shrink.seg.ai(p.arm)
					p.arm[nrow(p.arm),4] <- p.max
				}
			}	
			if(nrow(p.arm) >= 2){
				p.arm <- cbind(p.arm[,1:8], c(0,1)[match((p.arm[,4]-p.arm[,3]) >= 10e6, c('FALSE','TRUE'))])
				for(k in 2:nrow(p.arm)){
					if(p.arm[k,9] == 1 & p.arm[(k-1),9]==1){
						sample.lst <- c(sample.lst, 1)
						if(return.loc){
							a<- cbind(p.arm[(k-1),1:8], 2,'p-arm') ## Number indicates if start (1) or end (2) defines the breakpoint
							b <- cbind(p.arm[k,1:8], 1,'p-arm')
							colnames(a)[9:10]<- colnames(b)[9:10]<- c('LST breakpoint', 'chr. arm')
							out.seg <- rbind(out.seg, a,b)
						}
					}
				}
			}
			if(nrow(q.arm) >= 2){
				q.arm <- cbind(q.arm[,1:8], c(0,1)[match((q.arm[,4]-q.arm[,3]) >= 10e6, c('FALSE','TRUE'))])
				for(k in 2:nrow(q.arm)){
					if(q.arm[k,9] == 1 & q.arm[(k-1),9]==1){
						sample.lst <- c(sample.lst, 1)
						if(return.loc){
							a<- cbind(q.arm[(k-1),1:8], 2,'q-arm') ## Number indicates if start (1) or end (2) defines the breakpoint
							b <- cbind(q.arm[k,1:8], 1,'q-arm')
							colnames(a)[9:10]<- colnames(b)[9:10]<- c('LST breakpoint', 'chr. arm')
							out.seg <- rbind(out.seg, a,b)
						}
					}
				}
			}
		}
		output[j] <- sum(sample.lst)
	}
	if(return.loc){
		return(out.seg)
	} else {
		return(output)
	}
}
		

##################
#This function is an implementation of the cisplatin predictor developed by Myriad with Gordon Mills PMID: 23047548.
#Abkevichs cutoffs: > 10, found in the supplementary info
##################
#' HRD-LOH index
#' 
#' This function is an implementation of the cisplatin predictor developed by Myriad with Gordon Mills PMID: 23047548.
#' Abkevichs cutoffs: > 10, found in the supplementary info.
#' 
#' @param seg seg must be an ASCAT output object, in DNAcopy format. Requires the following columns in the given order: "SampleID", "Chromosome", "Start", "End", "nProbes", "totalCN", "nA", "nB", "Ploidy" and	"AberrantCellFraction". Column "nProbes", "Ploidy" and	"AberrantCellFraction" are not used here and can be set to NA.
#' @param nA is the column where copy number of A allele is found. This needs to be set to 7, because of the columns are referenced by position according to the required column order of seg.
#' @author  Nicolai Juul Birkbak, njuul@cbs.dtu.dk
#' @references Abkevich, V., Timms, K. M., Hennessy, B. T., Potter, J., Carey, M. S., Meyer, L. a., ... Lanchbury, J. S. (2012). Patterns of genomic loss of heterozygosity predict homologous recombination repair defects in epithelial ovarian cancer. British Journal of Cancer, 107(10), 1776–82. https://doi.org/10.1038/bjc.2012.451
#' @export
#' @examples 
#' HRD_LOH <- calc.hrd(ascat.data2, nA=7,check.names=FALSE, return.loc=FALSE)
calc.hrd <- function(seg, nA=7,check.names=FALSE, return.loc=FALSE){
	#Seg must be an ASCAT output object, in DNAcopy format. 
	#nA is the column where copy number of A allele is found
	#edit 20140225: added a check names function
	#edit 20140226: added return.loc to return breakpoint borders.
	seg[,1] <- as.character(seg[,1])
	if(check.names){
		seg <- check.names.fn(seg)
	}
	if(! all(seg[,8] <= seg[,7]) ){
		cat("Warning!! nB  not always <= nA!!  -- Correcting for internal use (only!)\n") # In case ASCAT people change the algorithm  ### They did!!
		tmp <- seg
		seg[tmp[,8] > tmp[,7],7]  <- tmp[tmp[,8] > tmp[,7],8]
		seg[tmp[,8] > tmp[,7],8]  <- tmp[tmp[,8] > tmp[,7],7]
	}
	nB <- nA+1
	output <- rep(0, length(unique(seg[,1])))
	names(output) <- unique(seg[,1])
	if(return.loc) {
			out.seg <- matrix(0,0,9)
			colnames(out.seg) <- c(colnames(seg)[1:8],'HRD breakpoint')
	}
	for(i in unique(seg[,1])){
		segSamp <- seg[seg[,1] %in% i,]
		chrDel <-vector() 
		for(j in unique(segSamp[,2])){ 	
			if(all(segSamp[segSamp[,2] == j,nB] == 0)) {
				chrDel <- c(chrDel, j)
			}
		}
		segLOH <- segSamp[segSamp[,nB] == 0 & segSamp[,nA] != 0,,drop=F]
		segLOH <- segLOH[segLOH[,4]-segLOH[,3] > 15e6,,drop=F]
		segLOH <- segLOH[!segLOH[,2] %in% chrDel,,drop=F]
		output[i] <- nrow(segLOH)
		if(return.loc){
			if(nrow(segLOH) < 1){next}
			segLOH <- cbind(segLOH[,1:8], 1)
			colnames(segLOH)[9] <- 'HRD breakpoint'
			out.seg <- rbind(out.seg, segLOH)
		}
	}
	if(return.loc){
		return(out.seg)
	} else {
		return(output)
	}
}	




################
#Number of telomeric AI, v1. 
#This is the original function as published. It has issues when samples are triploid. I suggest using the new version, calc.ai
################
#Edit 29/4-2012: To remove dependency on chrominfo, in the case there is no such matrix supplied (relevant when no centromere is mapped, as in chickens).
no.tel <- function(seg.out, chrominfo=chrominfo, min.size=0, min.probes=500, max.size=1e9, cnv.check="no", cnv.seg=NULL, cnv.gain=NULL,check.names=FALSE){

	# seg.out = segmented output in the form of a 6 column matrix like the output of DNA copy. Sixth column being AI yes/no (1/0). This is done simply by taking the nA and nB column from ASCAT output, and asking if they are the same, turn true/false into 1/0, and put this into the 6th column.
	# chrominfo = a 3 column matrix with information about the chromosomes: chromosome name, chromosome length, centromere location.
	# min.size = minimum size of segments
	# min.probes = minimum number of probes in segments (I use 500 for the 900,000 probe SNP6, then scale down)
	# max.size = maximum segment size (usually this is not an issue)
	# cnv.check = check for copy number status of AI segments
	# cnv.seg = a segmented copy number matrix, DNAcopy style. Ignored unless cnv.check is not "no"
	# cnv.gain = cnv threshold larger than which AI should not be counted

	#edit 2013-04-15: added check.names
	if(check.names){
		seg <- check.names.fn(seg)
	}

	if(class(seg.out) == "DNAcopy"){
		seg.out <- seg.out$output
	}
	tmp.segs <- seg.out[!seg.out[,2] %in% c("MT", "Y", "24"),]
	tmp.segs <- tmp.segs[!tmp.segs[,5] < min.probes,]
	tmp.segs[,2] <- as.character(tmp.segs[,2])
	if(nrow(tmp.segs[tmp.segs[,2] == "X",]) > 0){
		tmp.segs[tmp.segs[,2] == "X",2] <- rep(23, nrow(tmp.segs[tmp.segs[,2] == "X",]))
	}
	if(cnv.check != "no"){
		if(class(cnv.seg) == "DNAcopy"){
			cnv.seg <- cnv.seg$output
		}
		tmp.cnv <- cnv.seg[!cnv.seg[,2] %in% c("MT", "Y", "24"),]
		tmp.cnv <- tmp.cnv[!tmp.cnv[,5] < min.probes,]
		tmp.cnv[,2] <- as.character(tmp.cnv[,2])
		if(nrow(tmp.cnv[tmp.cnv[,2] == "X",]) > 0){
			tmp.cnv[tmp.cnv[,2] == "X",2] <- rep(23, nrow(tmp.cnv[tmp.cnv[,2] == "X",]))
		}
	}
	if(!all(seg.out[,6] %in% c(0,1))){ stop("Segmented input must be in 0/1 format!")}
	tmp.segs[tmp.segs[,6] == 1,6] <- 2
	for(j in 1:length(unique(tmp.segs[,1]))){
		tmp.sample <- tmp.segs[tmp.segs[,1]== unique(tmp.segs[,1])[j],]
		for(i in unique(tmp.segs[,2])){
			tmp1 <- tmp.sample[tmp.sample[,2]==i,,drop=F]
			if(nrow(tmp1) == 0){ next}
			if(class(chrominfo) == 'logical'){			# By logical, we assume that chrominfo == FALSE, hence we assume it is not there.	
				if(tmp1[1,6] == 2 & nrow(tmp1) != 1){
					tmp.sample[tmp.sample[,2]==i,6][1] <- 1
					if(cnv.check != "no"){
						if(cnv.seg[cnv.seg[,1] == unique(tmp.segs[,1])[j] & cnv.seg[,2] == i,][1,6] > cnv.gain){
							tmp.sample[tmp.sample[,2]==i,6][1] <- 0
						}
					}
				}
				if(tmp1[nrow(tmp1),6] == 2 & nrow(tmp1) != 1){
					tmp.sample[tmp.sample[,2]==i,6][nrow(tmp.sample[tmp.sample[,2]==i,])] <- 1
					if(cnv.check != "no"){
						if(cnv.seg[cnv.seg[,1] == unique(tmp.segs[,1])[j] & cnv.seg[,2] == i,][nrow(cnv.seg[cnv.seg[,1] == unique(tmp.segs[,1])[j] & cnv.seg[,2] == i,]),6] > cnv.gain){
							tmp.sample[tmp.sample[,2]==i,6][nrow(tmp.sample[tmp.sample[,2]==i,])] <- 0
						}
					}
				}
			}
			if(class(chrominfo) != 'logical'){				
				if(tmp1[1,6] == 2 & nrow(tmp1) != 1 & tmp1[1,4] < (chrominfo[i,3]*1000)){
					tmp.sample[tmp.sample[,2]==i,6][1] <- 1
					if(cnv.check != "no"){
						if(cnv.seg[cnv.seg[,1] == unique(tmp.segs[,1])[j] & cnv.seg[,2] == i,][1,6] > cnv.gain){
							tmp.sample[tmp.sample[,2]==i,6][1] <- 0
						}
					}
				}
				if(tmp1[nrow(tmp1),6] == 2 & nrow(tmp1) != 1 & tmp1[nrow(tmp1),3] > (chrominfo[i,3]*1000)){
					tmp.sample[tmp.sample[,2]==i,6][nrow(tmp.sample[tmp.sample[,2]==i,])] <- 1
					if(cnv.check != "no"){
						if(cnv.seg[cnv.seg[,1] == unique(tmp.segs[,1])[j] & cnv.seg[,2] == i,][nrow(cnv.seg[cnv.seg[,1] == unique(tmp.segs[,1])[j] & cnv.seg[,2] == i,]),6] > cnv.gain){
							tmp.sample[tmp.sample[,2]==i,6][nrow(tmp.sample[tmp.sample[,2]==i,])] <- 0
						}
					}
				}
			}
			if(nrow(tmp.sample[tmp.sample[,2]==i,]) == 1 & tmp.sample[tmp.sample[,2]==i,6][1] != 0){
				tmp.sample[tmp.sample[,2]==i,6][1] <- 3
			}
		}
		tmp.segs[tmp.segs[,1]== unique(tmp.segs[,1])[j],] <- tmp.sample
	}	
	#0 = no LOH, 1=telomeric LOH, 2=interstitial LOH, 3= whole chromosome LOH
	no.events <- matrix(0, nrow=length(unique(tmp.segs[,1])), ncol=(5+length(unique(tmp.segs[,2]))))
	rownames(no.events) <- unique(tmp.segs[,1])
	colnames(no.events) <- c("Telomeric AI", "Mean size", "Interstitial AI", "Mean Size", "Whole chr AI", sort(as.numeric(unique(tmp.segs[,2]))))
	for(i in unique(tmp.segs[,1])){
		tmp <- tmp.segs[tmp.segs[,1] == i,]
		tmp <- tmp[(tmp[,4] - tmp[,3]) > min.size,]
		tmp <- tmp[(tmp[,4] - tmp[,3]) < max.size,]
		no.events[i,1] <- nrow(tmp[tmp[,6] == 1,])
		no.events[i,2] <- mean(tmp[tmp[,6] == 1,4] - tmp[tmp[,6] == 1,3])
		no.events[i,3] <- nrow(tmp[tmp[,6] == 2,])
		no.events[i,4] <- mean(tmp[tmp[,6] == 2,4] - tmp[tmp[,6] == 2,3])
		no.events[i,5] <- nrow(tmp[tmp[,6] == 3,])
		no.events[i,tmp[tmp[,6] == 3,2]] <- 1
	}
	return(no.events)
}



################
#Function to check for duplicate sample names in an ASCAT segmented object, and return new names
################
check.names.fn <- function(seg, max.size=3.2e9, remove.dup =TRUE){
	tmp <- setNames(paste(seg[,1],seg[,9],seg[,11],sep='_'), seg[,1])
	sample.names <-	names(tmp[!duplicated(tmp)])# this is based on the fact that ASCAT ploidy is always highly variable, as it is a mean of all the segments
	if(any(duplicated(sample.names))){
		cat("Warning!! Duplicate names found! Renaming duplicates, adding a '.(number)' to each beyond no. 1\n")
		dup.samp <- unique(sample.names[duplicated(sample.names)])
		for(i in 1:length(dup.samp)){
			dup.samp.seg <- seg[seg[,1] %in% dup.samp[i],]		
			tmp<- unique(paste(dup.samp.seg[,1],dup.samp.seg[,9],dup.samp.seg[,11],sep='_'))
			for(k in 2:length(tmp)){
				dup.samp.seg[paste(dup.samp.seg[,1],dup.samp.seg[,9],dup.samp.seg[,11],sep='_') %in% tmp[k],1] <- paste(dup.samp.seg[paste(dup.samp.seg[,1],dup.samp.seg[,9],dup.samp.seg[,11],sep='_') %in% tmp[k],1],'.',k,sep='')
			}
			seg[seg[,1] %in% dup.samp[i],] <- dup.samp.seg
		}
	}
	# Now check for "silent" duplicates
	sample.names <- unique(seg[,1])
	for(i in sample.names){
		tmp <- seg[seg[,1] %in% i,]
		if(sum(tmp[,4] - tmp[,3]) > max.size){
			a <- which(tmp[,2] %in% 1)
			b <- which(tmp[,2] %in% 22)
			a <- a[!(a == (1:b[1])[1:length(a)])][1]
			tmp[a:nrow(tmp),1] <- paste('Duplicate', tmp[1,1], sep='-')
			seg[seg[,1] %in% i,] <- tmp
		}
	}
	if(remove.dup){
		seg <- seg[!grepl('Duplicate', seg[,1]),]
	}
	return(seg)
}


################
#Function to shrink a DNAcopy style matrix. Works on a chromosome level. Used to condense rows with identical CN values, which may have been separated due to other values not currently of interest (e.g. allelic copy differences)
#Required by min.cons.new function above.
################
shrink.seg <- function(chr.seg)	{	
	new.chr <- matrix(0,0,6)
	for(i in unique(chr.seg[,1])){
		tmp <- chr.seg[chr.seg[,1] %in% i,]
		if(any(duplicated(tmp[,6]))){
			seg.class <- c(1)
			for(j in 2:nrow(tmp)){
				ifelse(tmp[(j-1),6] == tmp[j,6], seg.class <- c(seg.class, seg.class[j-1]), seg.class <- c(seg.class, seg.class[j-1]+1))
			}
			for(j in unique(seg.class)){
				tmp[seg.class %in% j,4] <- max(tmp[seg.class %in% j,4])
			}
			tmp<- tmp[!duplicated(seg.class),]
		}
		new.chr <- rbind(new.chr, tmp)
	}
	colnames(new.chr) <- colnames(chr.seg)
	return(new.chr)
}
		
		
################
#AI-specific version of above function to shrink a DNAcopy style matrix. Works on a chromosome level. Used to condense rows with identical CN values, which may have been separated due to other values not currently of interest (e.g. filtered out for number of probes)
################
shrink.seg.ai <- function(chr.seg)	{	
	new.chr <- matrix(0,0,ncol(chr.seg))
	colnames(new.chr) <- colnames(chr.seg)
	new.chr <- chr.seg
	seg.class <- c(1)
	for(j in 2:nrow(new.chr)){
		ifelse(new.chr[(j-1),7] == new.chr[j,7] & new.chr[(j-1),8] == new.chr[j,8], seg.class <- c(seg.class, seg.class[j-1]), seg.class <- c(seg.class, seg.class[j-1]+1))
	}
	for(j in unique(seg.class)){
		new.chr[seg.class %in% j,4] <- max(new.chr[seg.class %in% j,4])
		new.chr[seg.class %in% j,5] <- sum(new.chr[seg.class %in% j,5])
	}
	new.chr<- new.chr[!duplicated(seg.class),]
	return(new.chr)
}
