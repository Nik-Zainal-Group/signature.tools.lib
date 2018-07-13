#' HRDetect Pipeline
#' 
#' Run the HRDetect pipeline
#' 
#' @export
#' 
HRDetect_pipeline <- function(data_matrix,
                              SNV_vcf_files=NULL,
                              SNV_tab_files=NULL,
                              SNV_catalogues=NULL,
                              Indels_vcf_files=NULL,
                              CNV_vcf_files=NULL,
                              SV_vcf_files=NULL,
                              SV_catalogues=NULL,
                              nparallel=1){
  #if multiple parallel cores are used, set it here
  if(nparallel>1){
    doMC::registerDoMC(nparallel)
  }
  
  #check that the matrix has correct features (columns)
  col_needed <- c("del.mh.prop", "e.3", "SV3", "SV5", "hrd", "e.8")
  if (!length(intersect(col_needed,colnames(data_matrix)))==length(col_needed)){
    stop("incorrect data matrix columns specified, you need the following columns: \"del.mh.prop\", \"e.3\", \"SV3\", \"SV5\", \"hrd\", \"e.8\"")
  }
  
  samples_list <- rownames(data_matrix)
  
  #check whether SNV related columns are NA or incomplete and if so check whether the catalogues
  #are available for sig fit, and if not check whether the vcf files are available for building catalogues
  
  need_to_compute_exposures <- any(is.na(data_matrix[,c("e.3","e.8")]))
  
  if(need_to_compute_exposures){
    #find out which samples that have no exposures have vcf file or catalogue
    incomplete_samples_pos <- which(apply(data_matrix[,c("e.3","e.8")],1,function(x) any(is.na(x))))
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
    #now, check that if a sample has both catalogue and vcf file, there is no need to compute the catalogue
    incomplete_samples_with_vcfSNV <- setdiff(incomplete_samples_with_vcfSNV,incomplete_samples_with_catalogueSNV)
    
    #now compute the 
  }
  
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
applyHRDetectDavies2017 <- function(data_matrix,features_names=c("del.mh.prop","e.3","SV3","SV5","hrd","e.8"),attachContributions=FALSE){
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
  
  hrdetect_input <- data_matrix[,features_names]
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
