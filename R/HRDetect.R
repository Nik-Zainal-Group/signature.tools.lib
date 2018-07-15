#' HRDetect Pipeline
#' 
#' Run the HRDetect pipeline
#' 
#' @export
#' 
HRDetect_pipeline <- function(data_matrix,
                              genome.v="hg19",
                              SNV_vcf_files=NULL,
                              SNV_tab_files=NULL,
                              SNV_catalogues=NULL,
                              Indels_vcf_files=NULL,
                              CNV_tab_files=NULL,
                              SV_bedpe_files=NULL,
                              SV_catalogues=NULL,
                              nparallel=1){
  #if multiple parallel cores are used, set it here
  if(nparallel>1){
    doMC::registerDoMC(nparallel)
  }
  
  #check that the matrix has correct features (columns)
  col_needed <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
  if (!length(intersect(col_needed,colnames(data_matrix)))==length(col_needed)){
    stop("[error HRDetect_pipeline] incorrect data_matrix columns specified, you need the following columns: \"del.mh.prop\", \"SNV3\", \"SV3\", \"SV5\", \"hrd\", \"SNV8\"")
  }
  # else{
  #   message("input data_matrix is formatted correctly")
  # }
  
  samples_list <- rownames(data_matrix)
  
  #check whether SNV related columns are NA or incomplete and if so check whether the catalogues
  #are available for sig fit, and if not check whether the vcf (or tab) files are available for building catalogues
  
  message("[info HRDetect_pipeline] Single Nucleotide Variations")
  
  SNV_cols <- c("SNV3","SNV8")
  
  need_to_compute_SNVexposures <- any(is.na(data_matrix[,SNV_cols]))
  
  if(need_to_compute_SNVexposures){
    
    message("[info HRDetect_pipeline] Some samples in the input data_matrix do not have the exposures for SNV3 and SNV8, checking if the user supplied SNV catalogues, VCF files or TAB files for those samples.")
    
    #find out which samples that have no exposures have vcf file or catalogue
    incomplete_samples_pos <- which(apply(data_matrix[,SNV_cols,drop=FALSE],1,function(x) any(is.na(x))))
    incomplete_samples <- rownames(data_matrix)[incomplete_samples_pos]
    if (!is.null(SNV_catalogues)){
      incomplete_samples_with_catalogueSNV <- intersect(incomplete_samples,colnames(SNV_catalogues))
    }else{
      #there is no SNV catalogue given, so no incomplete sample has a catalogue
      incomplete_samples_with_catalogueSNV <- character(0)
    }
    if (!is.null(SNV_vcf_files)){
      incomplete_samples_with_vcfSNV <- intersect(incomplete_samples,names(SNV_vcf_files))
    }else{
      #there is no SNV vcf files given, so no incomplete sample has a vcf file
      incomplete_samples_with_vcfSNV <- character(0)
    }
    if (!is.null(SNV_tab_files)){
      incomplete_samples_with_tabSNV <- intersect(incomplete_samples,names(SNV_tab_files))
    }else{
      #there is no SNV tab files given, so no incomplete sample has a tab file
      incomplete_samples_with_tabSNV <- character(0)
    }
    #now, check that if a sample has both catalogue and vcf (or tab) file, there is no need to compute the catalogue
    incomplete_samples_with_vcfSNV <- setdiff(incomplete_samples_with_vcfSNV,incomplete_samples_with_catalogueSNV)
    incomplete_samples_with_tabSNV <- setdiff(incomplete_samples_with_tabSNV,incomplete_samples_with_catalogueSNV)
    #also if a sample has both vcf and tab file, use the vcf file
    incomplete_samples_with_tabSNV <- setdiff(incomplete_samples_with_tabSNV,incomplete_samples_with_vcfSNV)
    
    #initialise the SNV catalogues data frame if necessary
    mut.order <- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")
    if(is.null(SNV_catalogues)){
      SNV_catalogues <- data.frame(row.names = mut.order)
    }else{
      SNV_catalogues <- SNV_catalogues[mut.order,,drop=FALSE]
    }
    
    #compute the SNV catalogue of samples with VCF files where necessary
    if (length(incomplete_samples_with_vcfSNV)>0){
      
      message("[info HRDetect_pipeline] VCF files will be converted to SNV catalogues for the following samples: ",paste(incomplete_samples_with_vcfSNV,collapse = " "))
      
      cat_list <- foreach::foreach(sample=incomplete_samples_with_vcfSNV) %dopar% {
        res <- vcfToSNVcatalogue(SNV_vcf_files[sample],genome.v = genome.v)
        colnames(res$catalogue) <- sample
        res$catalogue
      }
      #add new SNV catalogues to the catalogues matrix
      for (i in 1:length(cat_list)){
        newcat <- cat_list[[i]]
        SNV_catalogues <- cbind(SNV_catalogues,newcat)
        incomplete_samples_with_catalogueSNV <- c(incomplete_samples_with_catalogueSNV,colnames(newcat))
      }
    }
    
    #compute the SNV catalogue of samples with TAB files where necessary
    if (length(incomplete_samples_with_tabSNV)>0){
      
      message("[info HRDetect_pipeline] TAB files will be converted to SNV catalogues for the following samples: ",paste(incomplete_samples_with_tabSNV,collapse = " "))
      
      cat_list <- foreach::foreach(sample=incomplete_samples_with_tabSNV) %dopar% {
        subs <- read.table(file = SNV_tab_files[sample],
                  sep = "\t",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
        res <- tabToSNVcatalogue(subs,genome.v = genome.v)
        colnames(res$catalogue) <- sample
        res$catalogue
      }
      #add new SNV catalogues to the catalogues matrix
      for (i in 1:length(cat_list)){
        newcat <- cat_list[[i]]
        SNV_catalogues <- cbind(SNV_catalogues,newcat)
        incomplete_samples_with_catalogueSNV <- c(incomplete_samples_with_catalogueSNV,colnames(newcat))
      }
    }
    
    #compute exposures for the samples in the catalogueSNV that do not have exposures yet
    #consider only the catalogues needed
    SNV_catalogues_toFit <- SNV_catalogues[,incomplete_samples_with_catalogueSNV,drop=FALSE]
    
    #run sigfit with bootstrap if there are samples in the SNV_catalogue
    if(length(incomplete_samples_with_catalogueSNV)>0){
      message("[info HRDetect_pipeline] COSMIC30 signatures exposures will be estiamated for the following samples: ",paste(incomplete_samples_with_catalogueSNV,collapse = " "))
      message("[info HRDetect_pipeline] Running Signature fit with 100 bootstraps. Increase sparsity by removing exposures with 5% threshold of total mutations and 0.05 threshold of p-value, i.e. exposure of a signature in a sample is set to zero if the probability of having less than 5% of total mutations assigned to that signature is greather than 0.05.")
      res <- SignatureFit_withBootstrap(SNV_catalogues_toFit,
                                        cosmic30,
                                        nboot = 100,
                                        threshold_percent = 5,
                                        threshold_p.value = 0.05,
                                        verbose = FALSE,
                                        nparallel = nparallel)
      #add the resulting exposures to the 
      res$E_median_filtered
      data_matrix[colnames(res$E_median_filtered),SNV_cols] <- t(res$E_median_filtered[c("Signature.3","Signature.8"),])
    }
    
  }
  
  #check whether SV related columns are NA or incomplete and if so check whether the catalogues
  #are available for sig fit, and if not check whether the BEDPE files are available for building catalogues
  
  message("[info HRDetect_pipeline] Structural Variants (Rearrangements)")
  
  SV_cols <- c("SV3","SV5")
  
  need_to_compute_SVexposures <- any(is.na(data_matrix[,SV_cols]))
  
  
  if(need_to_compute_SVexposures){
    
    message("[info HRDetect_pipeline] Some samples in the input data_matrix do not have the exposures for SV3 and SV5, checking if the user supplied SV catalogues or BEDPE files for those samples.")
    
    #find out which samples that have no exposures have BEDPE file or catalogue
    incomplete_samples_pos <- which(apply(data_matrix[,SV_cols,drop=FALSE],1,function(x) any(is.na(x))))
    incomplete_samples <- rownames(data_matrix)[incomplete_samples_pos]
    if (!is.null(SV_catalogues)){
      incomplete_samples_with_catalogueSV <- intersect(incomplete_samples,colnames(SV_catalogues))
    }else{
      #there is no SV catalogue given, so no incomplete sample has a catalogue
      incomplete_samples_with_catalogueSV <- character(0)
    }
    if (!is.null(SV_bedpe_files)){
      incomplete_samples_with_bedpeSV <- intersect(incomplete_samples,names(SV_bedpe_files))
    }else{
      #there is no SV bedpe files given, so no incomplete sample has a bedpe file
      incomplete_samples_with_bedpeSV <- character(0)
    }
    
    #now, check that if a sample has both catalogue and BEDPE file, there is no need to compute the catalogue
    incomplete_samples_with_bedpeSV <- setdiff(incomplete_samples_with_bedpeSV,incomplete_samples_with_catalogueSV)
    
    #initialise the SV catalogues data frame if necessary
    catalogue.labels <- c('clustered_del_1-10Kb', 'clustered_del_10-100Kb', 'clustered_del_100Kb-1Mb', 'clustered_del_1Mb-10Mb', 'clustered_del_>10Mb', 'clustered_tds_1-10Kb', 'clustered_tds_10-100Kb', 'clustered_tds_100Kb-1Mb', 'clustered_tds_1Mb-10Mb', 'clustered_tds_>10Mb', 'clustered_inv_1-10Kb', 'clustered_inv_10-100Kb', 'clustered_inv_100Kb-1Mb', 'clustered_inv_1Mb-10Mb', 'clustered_inv_>10Mb', 'clustered_trans', 'non-clustered_del_1-10Kb', 'non-clustered_del_10-100Kb', 'non-clustered_del_100Kb-1Mb', 'non-clustered_del_1Mb-10Mb', 'non-clustered_del_>10Mb', 'non-clustered_tds_1-10Kb', 'non-clustered_tds_10-100Kb', 'non-clustered_tds_100Kb-1Mb', 'non-clustered_tds_1Mb-10Mb', 'non-clustered_tds_>10Mb', 'non-clustered_inv_1-10Kb', 'non-clustered_inv_10-100Kb', 'non-clustered_inv_100Kb-1Mb', 'non-clustered_inv_1Mb-10Mb', 'non-clustered_inv_>10Mb', 'non-clustered_trans')
    if(is.null(SV_catalogues)){
      SV_catalogues <- data.frame(row.names = catalogue.labels)
    }else{
      SV_catalogues <- SV_catalogues[catalogue.labels,,drop=FALSE]
    }
    
    #compute the SV catalogue of samples with BEDPE files where necessary
    if (length(incomplete_samples_with_bedpeSV)>0){
      
      message("[info HRDetect_pipeline] BEDPE files will be converted to Rearrangement catalogues for the following samples: ",paste(incomplete_samples_with_bedpeSV,collapse = " "))
      
      cat_list <- foreach::foreach(sample=incomplete_samples_with_bedpeSV) %dopar% {
        sv_bedpe <- read.table(SV_bedpe_files[sample],sep = "\t",header = TRUE,
                               stringsAsFactors = FALSE,check.names = FALSE)
        res <- bedpeToRearrCatalogue(sv_bedpe)
        colnames(res) <- sample
        res
      }
      #add new SV catalogues to the catalogues matrix
      for (i in 1:length(cat_list)){
        newcat <- cat_list[[i]]
        SV_catalogues <- cbind(SV_catalogues,newcat)
        incomplete_samples_with_catalogueSV <- c(incomplete_samples_with_catalogueSV,colnames(newcat))
      }
    }
    
    #compute exposures for the samples in the catalogueSV that do not have exposures yet
    #consider only the catalogues needed
    SV_catalogues_toFit <- SV_catalogues[,incomplete_samples_with_catalogueSV,drop=FALSE]
    
    #run sigfit with bootstrap if there are samples in the SV_catalogue
    if(length(incomplete_samples_with_catalogueSV)>0){
      message("[info HRDetect_pipeline] Breast 560 rearrangement signatures exposures will be estiamated for the following samples: ",paste(incomplete_samples_with_catalogueSV,collapse = " "))
      message("[info HRDetect_pipeline] Running Signature fit with 100 bootstraps. Increase sparsity by removing exposures with 5% threshold of total mutations and 0.05 threshold of p-value, i.e. exposure of a signature in a sample is set to zero if the probability of having less than 5% of total mutations assigned to that signature is greather than 0.05.")
      res <- SignatureFit_withBootstrap(SV_catalogues_toFit,
                                        RS.Breast560,
                                        nboot = 100,
                                        threshold_percent = 5,
                                        threshold_p.value = 0.05,
                                        verbose = FALSE,
                                        nparallel = nparallel)
      #add the resulting exposures to the 
      res$E_median_filtered
      data_matrix[colnames(res$E_median_filtered),SV_cols] <- t(res$E_median_filtered[c("RS3","RS5"),])
    }
    
  }
  
  #check whether small deletion at micro-homology (MH) related column (del.mh.prop) is NA or incomplete, and if so 
  #check whether the VCF indel files are available for computing the proportion of indels with MH
  
  message("[info HRDetect_pipeline] Deletions at Micro-homology (Indels)")
  
  MH_cols <- c("del.mh.prop")
  
  need_to_compute_MH <- any(is.na(data_matrix[,MH_cols]))
  
  
  if(need_to_compute_MH){
    
    message("[info HRDetect_pipeline] Some samples in the input data_matrix do not have the proportion of deletions at micro-homology, checking if the user supplied VCF indels files for those samples.")
    
    #find out which samples that have no del.mh.prop have vcf indel file
    incomplete_samples_pos <- which(apply(data_matrix[,MH_cols,drop=FALSE],1,function(x) any(is.na(x))))
    incomplete_samples <- rownames(data_matrix)[incomplete_samples_pos]
    if (!is.null(Indels_vcf_files)){
      incomplete_samples_with_vcfIndels <- intersect(incomplete_samples,names(Indels_vcf_files))
    }else{
      #there is no Indels vcf files given, so no incomplete sample has a vcf file
      incomplete_samples_with_vcfIndels <- character(0)
    }
    
    #compute del.mh.prop for the samples with vcf Indels file that do not have del.mh.prop yet
    if(length(incomplete_samples_with_vcfIndels)>0){
      message("[info HRDetect_pipeline] Proportion of Indels with MH will be computed for the following samples: ",paste(incomplete_samples_with_vcfIndels,collapse = " "))
      
      mh_list <- foreach::foreach(sample=incomplete_samples_with_vcfIndels) %dopar% {
        res <- vcfToIndelsClassification(Indels_vcf_files[sample],sample,genome.v = genome.v)
        new.prop.del.mh <- res$count_proportion["del.mh.prop"]
        new.prop.del.mh <- unlist(new.prop.del.mh)
        names(new.prop.del.mh) <- sample
        new.prop.del.mh
      }
      #add resulting del.mh.prop to the data_matrix
      for (i in 1:length(mh_list)){
        new.prop.del.mh <- mh_list[[i]]
        data_matrix[names(new.prop.del.mh),MH_cols] <- new.prop.del.mh
      }
    }
    
  }
  
  #check whether HRD-LOH related column (hrd) is NA or incomplete, and if so 
  #check whether the tab CNV files are available for computing the HRD-LOH index
  
  message("[info HRDetect_pipeline] HRD-LOH index (CNV)")
  
  hrd_cols <- c("hrd")
  
  need_to_compute_hrd <- any(is.na(data_matrix[,hrd_cols]))
  
  
  if(need_to_compute_hrd){
    
    message("[info HRDetect_pipeline] Some samples in the input data_matrix do not have the HRD-LOH index, checking if the user supplied tab CNV files for those samples.")
    
    #find out which samples that have no HRD-LOH have tab CNV file
    incomplete_samples_pos <- which(apply(data_matrix[,hrd_cols,drop=FALSE],1,function(x) any(is.na(x))))
    incomplete_samples <- rownames(data_matrix)[incomplete_samples_pos]
    if (!is.null(CNV_tab_files)){
      incomplete_samples_with_tabCNV <- intersect(incomplete_samples,names(CNV_tab_files))
    }else{
      #there is no CNV tab files given, so no incomplete sample has a CNV tab file
      incomplete_samples_with_tabCNV <- character(0)
    }
    
    #compute HRD-LOH for the samples with CNV tab file that do not have hrd yet
    if(length(incomplete_samples_with_tabCNV)>0){
      message("[info HRDetect_pipeline] HRD-LOH will be computed for the following samples: ",paste(incomplete_samples_with_tabCNV,collapse = " "))
      
      hrd_list <- foreach::foreach(sample=incomplete_samples_with_tabCNV) %dopar% {
        ascat.df <- read.table(CNV_tab_files[sample],sep = "\t",header = TRUE,
                               stringsAsFactors = FALSE,check.names = FALSE)
        res <- ascatToHRDLOH(ascat.df,sample)
        names(res) <- sample
        res
      }
      #add resulting hrd to the data_matrix
      for (i in 1:length(hrd_list)){
        res <- hrd_list[[i]]
        data_matrix[names(res),hrd_cols] <- res
      }
    }
    
  }
  
  #now, compute HRDetect score for all complete cases, so first notify of incomplete cases if any
  incomplete_cases <- rownames(data_matrix)[!complete.cases(data_matrix)]
  if(length(incomplete_cases)>0){
    message("[info HRDetect_pipeline] Some samples do not have data necessary to compute the HRDetect score (check output $data_matrix). Will not compute HRDetect score for the following samples: ",paste(incomplete_cases,collapse = " "))
    hrdetect_input <- data_matrix[complete.cases(data_matrix),,drop=FALSE]
  }else{
    hrdetect_input <- data_matrix
  }
  
  hrdetect_output <- applyHRDetectDavies2017(hrdetect_input, attachContributions = TRUE)
  
  
  #--- return results ---
  
  res <- list()
  res$data_matrix <- data_matrix
  res$SV_catalogues <- SV_catalogues
  res$SNV_catalogues <- SNV_catalogues
  return(res)
  
}



#' Apply HRDetect (Davies et al. 2017)
#' 
#' Apply HRDetect with the coeffcients from the Davies et al. 2017 publication.
#' This function requires a data frame with six features (columns) for each sample (rows),
#' and a list of features indicating which columns of the given matrix correspond to the required features.
#' Required features are: 
#' 1) proportion of deletions with microhomology, 
#' 2) number of mutations of substitution signature 3,
#' 3) number of mutations of rearrangemet signature 3,
#' 4) number of mutations of rearrangemet signature 5,
#' 5) HRD LOH index,
#' 6) number of mutations of substitution signature 8.
#' The function will return the HRDetect BRCAness probabilities, along with (optionally) the contributions
#' of each of the six features in each samples. The contributions are the normalised (log transoform and standardise) 
#' values of the features multiplied for the corresponding HRDetect logistic model coefficient.
#' 
#' @param data_matrix data frame containing a row for each sample and at least the columns specified in the features_names parameter
#' @param features_names list of column names of the matrix data_matrix. These indicate the features to be passed to HRDetect and should correspond to (in the exact order): 
#' 1) proportion of deletions with microhomology, 
#' 2) number of mutations of substitution signature 3,
#' 3) number of mutations of rearrangemet signature 3,
#' 4) number of mutations of rearrangemet signature 5,
#' 5) HRD LOH index,
#' 6) number of mutations of substitution signature 8.
#' @param attachContributions set to TRUE if you would like to have the contributions of the individual features to the samples HRDetect BRCAness probabilities.
#' @export
#' @references Davies, H., Glodzik, D., Morganella, S., Yates, L. R., Staaf, J., Zou, X., ... Nik-Zainal, S. (2017). HRDetect is a predictor of BRCA1 and BRCA2 deficiency based on mutational signatures. Nature Medicine, 23(4), 517–525. https://doi.org/10.1038/nm.4292
#' @examples
#' BRCAprob <- applyHRDetectDavies2017(data_matrix,features_names)
applyHRDetectDavies2017 <- function(data_matrix,features_names=c("del.mh.prop","SNV3","SV3","SV5","hrd","SNV8"),attachContributions=FALSE){
  #data_matrix:    rows are samples and columns are features
  #features_names: give the names of the features in data_matrix that correspond, in order, to:
  #                1) proportion of deletions with microhomology
  #                2) number of mutations of substitution signature 3
  #                3) number of mutations of rearrangemet signature 3
  #                4) number of mutations of rearrangemet signature 5
  #                5) HRD LOH index
  #                6) number of mutations of substitution signature 8
  
  
  features_mean <- c(0.218,
                     2.096,
                     1.260,
                     1.935,
                     2.195,
                     4.390)
  
  features_sd <- c(0.090,
                   3.555,
                   1.657,
                   1.483,
                   0.750,
                   3.179)
  
  features_weight <- c(2.398,
                       1.611,
                       1.153,
                       0.847,
                       0.667,
                       0.091)
  
  intercept <- -3.364
  
  hrdetect_input <- data_matrix[,features_names,drop=FALSE]
  #log and scale
  hrdetect_input <- log(hrdetect_input + 1)
  hrdetect_input <- hrdetect_input - matrix(rep(features_mean,nrow(hrdetect_input)),nrow = nrow(hrdetect_input),byrow = TRUE)
  hrdetect_input <- hrdetect_input/matrix(rep(features_sd,nrow(hrdetect_input)),nrow = nrow(hrdetect_input),byrow = TRUE)
  
  #get the score!
  linear_part <- rep(intercept,nrow(hrdetect_input)) + apply(hrdetect_input*matrix(rep(features_weight,nrow(hrdetect_input)),nrow = nrow(hrdetect_input),byrow = TRUE),1,sum)
  BRCA_prob <- 1/(1+exp(-linear_part))
  BRCA_prob <- matrix(BRCA_prob,nrow = length(BRCA_prob),ncol = 1,dimnames = list(names(BRCA_prob),"Probability"))
  
  if (attachContributions){
    non_zero <- features_names
    contrib <- hrdetect_input[,non_zero,drop=FALSE] * matrix(features_weight,ncol = length(non_zero),nrow = nrow(hrdetect_input),byrow = TRUE)
    intercept_matrix <- matrix(intercept,ncol = 1,nrow = nrow(hrdetect_input),dimnames = list(rownames(hrdetect_input),"intercept"))
    BRCA_prob <- cbind(intercept_matrix,contrib,BRCA_prob)
  }
  
  return(BRCA_prob)
}
