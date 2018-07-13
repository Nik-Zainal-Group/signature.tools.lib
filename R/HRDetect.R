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
  col_needed <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
  if (!length(intersect(col_needed,colnames(data_matrix)))==length(col_needed)){
    stop("incorrect data matrix columns specified, you need the following columns: \"del.mh.prop\", \"SNV3\", \"SV3\", \"SV5\", \"hrd\", \"SNV8\"")
  }
  
  samples_list <- rownames(data_matrix)
  
  #check whether SNV related columns are NA or incomplete and if so check whether the catalogues
  #are available for sig fit, and if not check whether the vcf (or tab) files are available for building catalogues
  
  SNV_cols <- c("SNV3","SNV8")
  
  need_to_compute_exposures <- any(is.na(data_matrix[,SNV_cols]))
  
  if(need_to_compute_exposures){
    #find out which samples that have no exposures have vcf file or catalogue
    incomplete_samples_pos <- which(apply(data_matrix[,SNV_cols],1,function(x) any(is.na(x))))
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
    
    #initialise the SNV catalogues data frame if necessary
    mut.order <- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")
    if(is.null(SNV_catalogues)){
      SNV_catalogues <- data.frame(row.names = mut.order)
    }else{
      SNV_catalogues <- SNV_catalogues[mut.order,,drop=FALSE]
    }
    
    #compute the SNV catalogue of samples where necessary
    for (sample in incomplete_samples_with_vcfSNV){
      
    }
    
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
