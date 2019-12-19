#' HRDetect Pipeline
#' 
#' Run the HRDetect pipeline. This function allows for flexible input 
#' specification to the HRDetect pipeline that computes the HRDetect 
#' BRCAness probability score as published in Davies et al. 2017.
#' It requires an input data frame "data_matrix", which contains a sample
#' in each row and one of six necessary features in each column. 
#' The six features can be computed by the pipeline if the necessary input files are provided.
#' The six features are:
#' 1) proportion of deletions at microhomology (del.mh.prop), 
#' 2) number of mutations of substitution signature 3 (SNV3),
#' 3) number of mutations of rearrangemet signature 3 (SV3),
#' 4) number of mutations of rearrangemet signature 5 (SV5),
#' 5) HRD LOH index (hrd),
#' 6) number of mutations of substitution signature 8 (SNV8).
#' For example, if the HRD LOH index has already been calculated, these can be
#' added to the input data_matrix, or if the SNV catalogues have already been calculated,
#' these can be supplied using the SNV_catalogues parameter while setting SNV3 and SNV8 columns
#' as "NA". Also, it is possible to provide different data for different samples. For example,
#' one can provide SNV3 and SNV8 number of mutations for some samples in data_matrix, while setting
#' SNV3 and SNV8 to NA for other samples, and providing either SNV catalogues and/or SNV VCF 
#' files for these samples. The function will return the HRDetect BRCAness probability score for
#' all the samples for which enough data are available to calculate all six necessary features.
#' Along with the score, the contribution of each feature to the score will be provided. In addition,
#' an updated data_matrix and other other data that have been calculated during the execution of the pipeline
#' will be returned as well.
#' 
#' Single Nucleotide Variations. Columns in data_matrix relative to SNV are SNV3 and SNV8. Values
#' corresponding to number of SNV3 and SNV8 mutations in each sample can be provided in the data frame data_matrix.
#' Alternatively, an SNV_catalogue data frame can be used to provide 96-channel SNV catalogues for the samples
#' (96-channels as rows and samples as columns). The 30 consensus COSMIC signatures will be fitted to 
#' the catalogues using a bootstrapping approach (Huang et al. 2017) and estimates for SNV3 and SNV8 will be added to the data_matrix.
#' Alternatively, SNV_catalogues can be constructed providing a list of either SNV VCF files or SNV TAB files.
#' 
#' Structural Variants (Rearrangements). Columns in data_matrix relative to SV are SV3 and SV5. Values
#' corresponding to number of SV3 and SV5 rearrangements in each sample can be provided in the data frame data_matrix.
#' Alternatively, an SV_catalogue data frame can be used to provide 32-channel SV catalogues for the samples
#' (32-channels as rows and samples as columns). The 6 Breast Cancer Rearrangement signatures (Nik-Zainal et al. 2016) will be fitted to 
#' the catalogues using a bootstrapping approach (Huang et al. 2017) and estimates for SV3 and SV5 will be added to the data_matrix.
#' Alternatively, SV_catalogues can be constructed providing a list of SV BEDPE files.
#' 
#' Deletions at Micro-homology (Indels). The column in data_matrix corresponding to the proportion of deletions at micro-homology is del.mh.prop.
#' The proportion of deletions at micro-homology for the samples can be calculated by the pipeline if the user provides Indels VCF files.
#' 
#' HRD-LOH index (CNV). The column in data_matrix corresponding to the HRD-LOH index is hrd.
#' The HRD-LOH index for the samples can be calculated by the pipeline if the user provides copy numbers TAB files.
#' 
#' @param data_matrix data frame containing a sample for each row and the six necessary features as columns. Columns should be labelled with the following names: del.mh.prop, SNV3, SV3, SV5, hrd, SNV8. Row names of the data frame should correspond to the sample names. If the values of the features need to be computed, set them to NA and provide additional data (e.g. catalogues, VCF/BEDPE/TAB files as specified in this documentation page).
#' @param genome.v genome version to use when constructing the SNV catalogue and classifying indels. Set it to either "hg19" or "hg38".
#' @param SNV_vcf_files list of file names corresponding to SNV VCF files to be used to construct 96-channel substitution catalogues. This should be a named vector, where the names indicate the sample name, so that each file can be matched to the corresponding row in the data_matrix input. The files should only contain SNV and should already be filtered according to the user preference, as all SNV in the file will be used and no filter will be applied.
#' @param SNV_tab_files list of file names corresponding to SNV TAB files to be used to construct 96-channel substitution catalogues. This should be a named vector, where the names indicate the sample name, so that each file can be matched to the corresponding row in the data_matrix input. The files should only contain SNV and should already be filtered according to the user preference, as all SNV in the file will be used and no filter will be applied. The files should contain a header in the first line with the following columns: chr, position, REF, ALT.
#' @param SNV_catalogues data frame containing 96-channel substitution catalogues. A sample for each column and the 96-channels as rows. Row names should have the correct channel names (see for example tests/testthat/test.snv.tab) and the column names should be the sample names so that each catalogue can be matched with the corresponding row in the data_matrix input.
#' @param Indels_vcf_files list of file names corresponding to Indels VCF files to be used to classify Indels and compute the proportion of indels at micro-homology. This should be a named vector, where the names indicate the sample name, so that each file can be matched to the corresponding row in the data_matrix input. The files should only contain indels (no SNV) and should already be filtered according to the user preference, as all indels in the file will be used and no filter will be applied.
#' @param CNV_tab_files list of file names corresponding to CNV TAB files (similar to ASCAT format) to be used to compute the HRD-LOH index. This should be a named vector, where the names indicate the sample name, so that each file can be matched to the corresponding row in the data_matrix input. The files should contain a header in the first line with the following columns: 'seg_no', 'Chromosome', 'chromStart', 'chromEnd', 'total.copy.number.inNormal', 'minor.copy.number.inNormal', 'total.copy.number.inTumour', 'minor.copy.number.inTumour'
#' @param SV_bedpe_files list of file names corresponding to SV (Rearrangements) BEDPE files to be used to construct 32-channel rearrangement catalogues. This should be a named vector, where the names indicate the sample name, so that each file can be matched to the corresponding row in the data_matrix input. The files should contain a rearrangement for each row (two breakpoint positions should be on one row as determined by a pair of mates of paired-end sequencing) and should already be filtered according to the user preference, as all rearrangements in the file will be used and no filter will be applied. The files should contain a header in the first line with the following columns: "chrom1", "start1", "end1", "chrom2", "start2", "end2" and "sample" (sample name). In addition, either two columns indicating the strands of the mates, "strand1" (+ or -) and "strand2" (+ or -), or one column indicating the structural variant class, "svclass": translocation, inversion, deletion, tandem-duplication. The column "svclass" should correspond to (Sanger BRASS convention): inversion (strands +/- or -/+ and mates on the same chromosome), deletion (strands +/+ and mates on the same chromosome), tandem-duplication (strands -/- and mates on the same chromosome), translocation (mates are on different chromosomes)..
#' @param SV_catalogues data frame containing 32-channel substitution catalogues. A sample for each column and the 32-channels as rows. Row names should have the correct channel names (see for example tests/testthat/test.cat) and the column names should be the sample names so that each catalogue can be matched with the corresponding row in the data_matrix input.
#' @param signature_type either "COSMIC" or one of the following organs: "Biliary", "Bladder", "Bone_SoftTissue", "Breast", "Cervix", "CNS", "Colorectal", "Esophagus", "Head_neck", "Kidney", "Liver", "Lung", "Lymphoid", "Ovary", "Pancreas", "Prostate", "Skin", "Stomach", "Uterus"
#' @param bootstrap_scores perform HRDetect score with bootstrap. This requires mutations or catalogues for subs/rearr to compute the bootstrap fit, and indels mutations to bootstrap the indels classification. HRD-LOH can still be provided using the input data_matrix.
#' @param nparallel how many parallel threads to use.
#' @return return a list that contains $data_matrix (updated input data_matrix with additional computed features), $hrdetect_output (data frame with HRDetect BRCAness Probability and contribution of the features), $SNV_catalogues (input SNV_catalogues updated with additional computed substitution catalogues if any), $SV_catalogues (input SV_catalogues updated with additional computed rearrangement catalogues if any)
#' @export
#' @references Davies, H., Glodzik, D., Morganella, S., Yates, L. R., Staaf, J., Zou, X., ... Nik-Zainal, S. (2017). HRDetect is a predictor of BRCA1 and BRCA2 deficiency based on mutational signatures. Nature Medicine, 23(4), 517–525. https://doi.org/10.1038/nm.4292
#' @references Nik-Zainal, S., Davies, H., Staaf, J., Ramakrishna, M., Glodzik, D., Zou, X., ... Stratton, M. R. (2016). Landscape of somatic mutations in 560 breast cancer whole-genome sequences. Nature, 534(7605), 1–20. https://doi.org/10.1038/nature17676
#' @references Huang, X., Wojtowicz, D., & Przytycka, T. M. (2017). Detecting Presence Of Mutational Signatures In Cancer With Confidence. bioRxiv, (October). https://doi.org/10.1101/132597
#' 
HRDetect_pipeline <- function(data_matrix=NULL,
                              genome.v="hg19",
                              SNV_vcf_files=NULL,
                              SNV_tab_files=NULL,
                              SNV_catalogues=NULL,
                              Indels_vcf_files=NULL,
                              CNV_tab_files=NULL,
                              SV_bedpe_files=NULL,
                              SV_catalogues=NULL,
                              signature_type="COSMIC",
                              cosmic_siglist=NULL,
                              bootstrap_scores=FALSE,
                              nparallel=1){
  #if multiple parallel cores are used, set it here
  doMC::registerDoMC(nparallel)
  
  #check that the matrix has correct features (columns)
  col_needed <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
  
  if (!is.null(data_matrix)){
    if (!length(intersect(col_needed,colnames(data_matrix)))==length(col_needed)){
      stop("[error HRDetect_pipeline] incorrect data_matrix columns specified, you need the following columns: \"del.mh.prop\", \"SNV3\", \"SV3\", \"SV5\", \"hrd\", \"SNV8\"")
    }
  }else{
    #there is no data_matrix, just build it from files
    #get the sample names from all the file lists
    sample_names <- unique(unlist(sapply(list(SNV_vcf_files,
                                SNV_tab_files,
                                SNV_catalogues,
                                Indels_vcf_files,
                                CNV_tab_files,
                                SV_bedpe_files,
                                SV_catalogues), function(x) names(x))))
    if(!is.null(sample_names)){
      data_matrix <- matrix(NA,ncol = 6,nrow = length(sample_names),dimnames = list(sample_names,col_needed))
    }else{
      stop("[error HRDetect_pipeline] no input data provided.")
    }
  }
  samples_list <- rownames(data_matrix)
  
  #check that the signature_type option is valid
  #get available organs
  available_organs <- intersect(unique(sapply(strsplit(colnames(all_organ_sigs_subs),split="_"),function(x) paste(x[1:(length(x)-1)],collapse = "_"))),
  unique(sapply(strsplit(colnames(all_organ_sigs_rearr),split="_"),function(x) paste(x[1:(length(x)-1)],collapse = "_"))))
  
  if(!(signature_type %in% c("COSMIC",available_organs))){
    stop(paste0("[error HRDetect_pipeline] invalid signature_type ",signature_type,"."))
  }
  
  #use either COSMIC sigs or tissue-specific sigs
  if(signature_type=="COSMIC"){
    if (is.null(cosmic_siglist)) cosmic_siglist <- 1:30
    sigstofit_subs <- cosmic30[,cosmic_siglist,drop=FALSE]
    sigstofit_rearr <- RS.Breast560
  }else{
    sigstofit_subs <- all_organ_sigs_subs[,colnames(all_organ_sigs_subs)[grep(pattern = paste0("^",signature_type),colnames(all_organ_sigs_subs))]]
    sigstofit_rearr <- all_organ_sigs_rearr[,colnames(all_organ_sigs_rearr)[grep(pattern = paste0("^",signature_type),colnames(all_organ_sigs_rearr))]]
  }
  
  #initialise bootstrap fit variables to NULL
  bootstrap_fit_subs <- NULL
  bootstrap_fit_rearr <- NULL
  #initialise final exposures output table
  exposures_subs <- NULL
  exposures_rearr <- NULL
  #initialise indels classification table
  indels_classification_table <- NULL
  #initialise hrdetect bootstrap tables
  hrdetect_bootstrap_table <- NULL
  q_5_50_95 <- NULL
  
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
      message("[info HRDetect_pipeline] Substitution signatures exposures will be estiamated for the following samples: ",paste(incomplete_samples_with_catalogueSNV,collapse = " "))
      message("[info HRDetect_pipeline] Running Signature fit with 100 bootstraps. Increase sparsity by removing exposures with 5% threshold of total mutations and 0.05 threshold of p-value, i.e. exposure of a signature in a sample is set to zero if the probability of having less than 5% of total mutations assigned to that signature is greather than 0.05.")
      bootstrap_fit_subs <- SignatureFit_withBootstrap(SNV_catalogues_toFit,
                                        sigstofit_subs,
                                        nboot = 100,
                                        threshold_percent = 5,
                                        threshold_p.value = 0.05,
                                        verbose = FALSE,
                                        nparallel = nparallel)
      #add the resulting exposures to the data_matrix
      exposures_subs <- bootstrap_fit_subs$E_median_filtered
      exposures_subs[is.nan(exposures_subs)] <- 0
      
      if(signature_type=="COSMIC"){
        if("Signature.3" %in% rownames(exposures_subs)){
          data_matrix[colnames(exposures_subs),"SNV3"] <- exposures_subs["Signature.3",]
        }else{
          data_matrix[colnames(exposures_subs),"SNV3"] <- 0
        }
        if("Signature.8" %in% rownames(exposures_subs)){
          data_matrix[colnames(exposures_subs),"SNV8"] <- exposures_subs["Signature.8",]
        }else{
          data_matrix[colnames(exposures_subs),"SNV8"] <- 0
        }
        # data_matrix[colnames(exposures_subs),SNV_cols] <- t(exposures_subs[c("Signature.3","Signature.8"),])
      }else{
        #convert to reference signatures
        exposures_subs <- t(conversion_matrix_subs[colnames(sigstofit_subs),]) %*% exposures_subs
        exposures_subs <- exposures_subs[apply(exposures_subs, 1,sum)>0,,drop=FALSE]
        #add to data_matrix if present
        if("Ref.Sig.3" %in% rownames(exposures_subs)){
          data_matrix[colnames(exposures_subs),"SNV3"] <- exposures_subs["Ref.Sig.3",]
        }else{
          data_matrix[colnames(exposures_subs),"SNV3"] <- 0
        }
        if("Ref.Sig.8" %in% rownames(exposures_subs)){
          data_matrix[colnames(exposures_subs),"SNV8"] <- exposures_subs["Ref.Sig.8",]
        }else{
          data_matrix[colnames(exposures_subs),"SNV8"] <- 0
        }
        
      }
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
                               stringsAsFactors = FALSE,check.names = FALSE,comment.char = "")
        reslist <- bedpeToRearrCatalogue(sv_bedpe)
        res <- reslist$rearr_catalogue
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
      message("[info HRDetect_pipeline] Rearrangement signatures exposures will be estiamated for the following samples: ",paste(incomplete_samples_with_catalogueSV,collapse = " "))
      message("[info HRDetect_pipeline] Running Signature fit with 100 bootstraps. Increase sparsity by removing exposures with 5% threshold of total mutations and 0.05 threshold of p-value, i.e. exposure of a signature in a sample is set to zero if the probability of having less than 5% of total mutations assigned to that signature is greather than 0.05.")
      bootstrap_fit_rearr <- SignatureFit_withBootstrap(SV_catalogues_toFit,
                                        sigstofit_rearr,
                                        nboot = 100,
                                        threshold_percent = 5,
                                        threshold_p.value = 0.05,
                                        verbose = FALSE,
                                        nparallel = nparallel)
      #add the resulting exposures to the data_matrix
      exposures_rearr <- bootstrap_fit_rearr$E_median_filtered
      exposures_rearr[is.nan(exposures_rearr)] <- 0
      
      if(signature_type=="COSMIC"){
        data_matrix[colnames(exposures_rearr),SV_cols] <- t(exposures_rearr[c("RS3","RS5"),])
      }else{
        #convert to reference signatures
        exposures_rearr <- t(conversion_matrix_rearr[colnames(sigstofit_rearr),]) %*% exposures_rearr
        exposures_rearr <- exposures_rearr[apply(exposures_rearr, 1,sum)>0,,drop=FALSE]
        #add to data_matrix if present
        if("RefSig R3" %in% rownames(exposures_rearr)){
          data_matrix[colnames(exposures_rearr),"SV3"] <- exposures_rearr["RefSig R3",]
        }else{
          data_matrix[colnames(exposures_rearr),"SV3"] <- 0
        }
        if("RefSig R5" %in% rownames(exposures_rearr)){
          data_matrix[colnames(exposures_rearr),"SV5"] <- exposures_rearr["RefSig R5",]
        }else{
          data_matrix[colnames(exposures_rearr),"SV5"] <- 0
        }
        if("RefSig R9" %in% rownames(exposures_rearr)){
          data_matrix[colnames(exposures_rearr),"SV5"] <- data_matrix[colnames(exposures_rearr),"SV5"] + exposures_rearr["RefSig R9",]
        }
      }
      
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
        res$count_proportion
      }
      #combine in one table and add to data_matrix
      indels_classification_table <- do.call(rbind,mh_list)
      rownames(indels_classification_table) <- indels_classification_table[,"sample"]
      data_matrix[rownames(indels_classification_table),MH_cols] <- indels_classification_table[,"del.mh.prop"]
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
  
  if(nrow(hrdetect_input)>0){
    message("[info HRDetect_pipeline] Computing HRDetect score and feature contributions for the following samples: ",paste(rownames(hrdetect_input),collapse = " "))
    hrdetect_output <- applyHRDetectDavies2017(hrdetect_input, attachContributions = TRUE)
    
    #attempt to run HRDetect with bootstrap if requested
    if(bootstrap_scores){
      message("[info HRDetect_pipeline] HRDetect boostrap scores requested!")
      
      #check whether we have all we need and for which samples
      bootstrap_samples <- c()
      
      #bootstrap exposures are required
      if(!is.null(bootstrap_fit_subs)) bootstrap_samples <- colnames(bootstrap_fit_subs$E_median_filtered)
      if(!is.null(bootstrap_fit_rearr)) bootstrap_samples <- intersect(bootstrap_samples,colnames(bootstrap_fit_rearr$E_median_filtered))
      
      #indels classification tables are required
      if(!is.null(indels_classification_table)) bootstrap_samples <- intersect(bootstrap_samples,rownames(indels_classification_table))
      
      #If samples are in the hrdetect_input table then they also have the HRD-LOH score
      bootstrap_samples <- intersect(bootstrap_samples,rownames(hrdetect_input))
      
      if(length(bootstrap_samples)>0){
        message("[info HRDetect_pipeline] Computing HRDetect boostrap for the following samples: ",paste(bootstrap_samples,collapse = ", "))
        
        #run hrdetect bootstrap scores for the bootstrap_samples
        nboots_hr <- 1000
        
        hrdetect_bootstrap_table <- list()
        for (j in 1:nboots_hr){
          tmp_data_matrix <- data_matrix[bootstrap_samples,,drop=FALSE]
          tmp_data_matrix[1:nrow(tmp_data_matrix),1:ncol(tmp_data_matrix)] <- NA
          
          #1) sample del.mh.prop
          deletions_muts <- t(indels_classification_table[,c("del.mh.count","del.rep.count","del.none.count")])[,bootstrap_samples,drop=FALSE]
          ndeletions <- apply(deletions_muts,2,sum)
          samples_deletions_muts <- generateRandMuts(deletions_muts)
          samples_ndeletions <- apply(samples_deletions_muts,2,sum)
          samples_deletions_prop <- t(samples_deletions_muts/matrix(rep(samples_ndeletions,3),nrow = nrow(samples_deletions_muts),ncol = ncol(samples_deletions_muts),byrow = TRUE))
          colnames(samples_deletions_prop) <- c("del.mh.prop","del.rep.prop","del.none.prop")
          tmp_data_matrix[bootstrap_samples,"del.mh.prop"] <- samples_deletions_prop[bootstrap_samples,"del.mh.prop"]
          
          #2) sample HRD-LOH index
          tmp_data_matrix[bootstrap_samples,"hrd"] <- rpois(length(bootstrap_samples),data_matrix[bootstrap_samples,"hrd"])
          
          #3) sample subs exposures
          subs_boots <- bootstrap_fit_subs$boot_list
          s1 <- sample(length(subs_boots),size = 1)
          current_subs <- t(subs_boots[[s1]])
          current_subs[is.na(current_subs)] <- 0
          #sparsity correction 5%
          sel <- t(apply(current_subs,1,function(x) x<sum(x)*0.05))
          current_subs[sel] <- 0
          if(signature_type=="COSMIC"){
            if("Signature.3" %in% colnames(current_subs)){
              tmp_data_matrix[bootstrap_samples,"SNV3"] <- current_subs[bootstrap_samples,"Signature.3"]
            }else{
              tmp_data_matrix[bootstrap_samples,"SNV3"] <- 0
            }
            if("Signature.8" %in% colnames(current_subs)){
              tmp_data_matrix[bootstrap_samples,"SNV8"] <- current_subs[bootstrap_samples,"Signature.8"]
            }else{
              tmp_data_matrix[bootstrap_samples,"SNV8"] <- 0
            }
          }else{
            #convert to ref sig
            exposures_RefSigs <- t(as.matrix(current_subs) %*% as.matrix(conversion_matrix_subs[colnames(current_subs),]))
            exposures_RefSigs <- exposures_RefSigs[apply(exposures_RefSigs, 1,sum)>0,,drop=FALSE]
            if("Ref.Sig.3" %in% rownames(exposures_RefSigs)){
              tmp_data_matrix[bootstrap_samples,"SNV3"] <- exposures_RefSigs["Ref.Sig.3",bootstrap_samples]
            }else{
              tmp_data_matrix[bootstrap_samples,"SNV3"] <- 0
            }
            if("Ref.Sig.8" %in% rownames(exposures_RefSigs)){
              tmp_data_matrix[bootstrap_samples,"SNV8"] <- exposures_RefSigs["Ref.Sig.8",bootstrap_samples]
            }else{
              tmp_data_matrix[bootstrap_samples,"SNV8"] <- 0
            }
          }
          #4) sample rearr exposures
          rearr_boots <- bootstrap_fit_rearr$boot_list
          s1 <- sample(length(rearr_boots),size = 1)
          current_rearr <- t(rearr_boots[[s1]])
          current_rearr[is.na(current_rearr)] <- 0
          #sparsity correction 5%
          sel <- t(apply(current_rearr,1,function(x) x<sum(x)*0.05))
          current_rearr[sel] <- 0
          if(signature_type=="COSMIC"){
            tmp_data_matrix[bootstrap_samples,"SV3"] <- current_rearr[bootstrap_samples,"RS3"]
            tmp_data_matrix[bootstrap_samples,"SV5"] <- current_rearr[bootstrap_samples,"RS5"]
          }else{
            #convert to ref sig
            exposures_RefSigs <- t(as.matrix(current_rearr) %*% as.matrix(conversion_matrix_rearr[colnames(current_rearr),]))
            exposures_RefSigs <- exposures_RefSigs[apply(exposures_RefSigs, 1,sum)>0,,drop=FALSE]
            if("RefSig R3" %in% rownames(exposures_RefSigs)){
              tmp_data_matrix[bootstrap_samples,"SV3"] <- exposures_RefSigs["RefSig R3",bootstrap_samples]
            }else{
              tmp_data_matrix[bootstrap_samples,"SV3"] <- 0
            }
            if("RefSig R5" %in% rownames(exposures_RefSigs)){
              tmp_data_matrix[bootstrap_samples,"SV5"] <- exposures_RefSigs["RefSig R5",bootstrap_samples]
            }else{
              tmp_data_matrix[bootstrap_samples,"SV5"] <- 0
            }
            if("RefSig R9" %in% rownames(exposures_RefSigs)){
              tmp_data_matrix[bootstrap_samples,"SV5"] <- tmp_data_matrix[bootstrap_samples,"SV5"] + exposures_RefSigs["RefSig R9",bootstrap_samples]
            }
          }
          
          hrdetect_bootstrap_table[[j]] <- t(applyHRDetectDavies2017(data_matrix = tmp_data_matrix,attachContributions = FALSE))
        }
        
        hrdetect_bootstrap_table <- do.call(rbind,hrdetect_bootstrap_table)
        rownames(hrdetect_bootstrap_table) <- 1:nrow(hrdetect_bootstrap_table)
        q_5_50_95 <- t(apply(hrdetect_bootstrap_table,2,function(x) quantile(x,c(0.05,0.5,0.95))))
        
        message("[info HRDetect_pipeline] HRDetect with bootstrap successful!")
        
      }else{
        message("[info HRDetect_pipeline] Not enough data to run HRDetect boostrap scores.")
      }
    }
    
    message("[info HRDetect_pipeline] HRDetect pipeline completed!")
  }else{
    message("[info HRDetect_pipeline] Impossible to compute any HRDetect score, as no sample seems to have all the necessary data. Check output $data_matrix to see what has been computed with the data provided.")
    hrdetect_output <- NULL
  }
  
  
  #--- return results ---
  
  res <- list()
  res$data_matrix <- data_matrix
  res$SV_catalogues <- SV_catalogues
  res$SNV_catalogues <- SNV_catalogues
  res$hrdetect_output <- hrdetect_output
  res$exposures_subs <- exposures_subs
  res$exposures_rearr <- exposures_rearr
  res$bootstrap_fit_subs <- bootstrap_fit_subs
  res$bootstrap_fit_rearr <- bootstrap_fit_rearr
  res$indels_classification_table <- indels_classification_table
  res$hrdetect_bootstrap_table <- hrdetect_bootstrap_table
  res$q_5_50_95 <- q_5_50_95
  return(res)
  
}



#' Apply HRDetect (Davies et al. 2017)
#' 
#' Apply HRDetect with the coeffcients from the Davies et al. 2017 publication.
#' This function requires a data frame with six features (columns) for each sample (rows),
#' and a list of features indicating which columns of the given matrix correspond to the required features.
#' Required features are: 
#' 1) proportion of deletions with microhomology (del.mh.prop), 
#' 2) number of mutations of substitution signature 3 (SNV3),
#' 3) number of mutations of rearrangemet signature 3 (SV3),
#' 4) number of mutations of rearrangemet signature 5 (SV5),
#' 5) HRD LOH index (hrd),
#' 6) number of mutations of substitution signature 8 (SNV8).
#' The function will return the HRDetect BRCAness probabilities, along with (optionally) the contributions
#' of each of the six features in each samples. The contributions are the normalised (log transoform and standardise) 
#' values of the features multiplied for the corresponding HRDetect logistic model coefficient.
#' 
#' @param data_matrix data frame containing a row for each sample and at least the columns specified in the features_names parameter
#' @param features_names list of column names of the matrix data_matrix. These indicate the features to be passed to HRDetect and should correspond to (in the exact order): 
#' 1) proportion of deletions with microhomology (del.mh.prop), 
#' 2) number of mutations of substitution signature 3 (SNV3),
#' 3) number of mutations of rearrangemet signature 3 (SV3),
#' 4) number of mutations of rearrangemet signature 5 (SV5),
#' 5) HRD LOH index (hrd),
#' 6) number of mutations of substitution signature 8 (SNV8).
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



#------------------------

#' plot_HRDLOH_HRDetect_Contributions
#' 
#' Function for plotting HRDetect results and contributions to HRDetect score.
#' 
#' @param file_name name of the output file (jpg)
#' @param HRDLOH_index HRD index score 
#' @param hrdetect_output output of HRDetect, containing contributions of each feature
#' @export
plot_HRDLOH_HRDetect_Contributions <- function(file_name,HRDLOH_index,hrdetect_output){
  
  col_needed <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
  
  #details plot
  hrdetect_output <- as.data.frame(hrdetect_output)
  reorder <- order(hrdetect_output$Probability,decreasing = TRUE)
  nsamples <- nrow(hrdetect_output)
  contributions <- hrdetect_output
  jpeg(filename = file_name,
       width = max(200+80*nsamples,2400),
       height = 1200,
       res = 200)
  par(mfrow=c(1,1))
  mat <- matrix(c(1,2,3),ncol = 1)
  layout(mat, c(1), c(0.8,0.65,1))
  
  par(mar=c(1,6,3,4))
  barplot(height = HRDLOH_index[reorder],
          ylim = c(0,max(HRDLOH_index)*1.1),
          ylab = "HRD-LOH score",
          cex.axis = 1.2,
          cex = 1.5,
          cex.lab=1.5)
  title(main=paste0("HRD-LOH score and BRCA1/BRCA2 deficiency score"), cex.main=2)
  abline(a=15,b=0,col="red")
  
  par(mar=c(1,6,1,4))
  barplot(height = hrdetect_output$Probability[reorder],
          ylim = c(0,1),ylab = "BRCA1/BRCA2\n def. score",
          names.arg = "",
          cex.axis = 1.2,
          cex.lab=1.5,
          cex.names = 1.5)
  abline(a=0.7,b=0,col="red")
  #text(x = -0.05,y = 0.75,labels = "0.7",col = "red")
  matrix_to_plot <- t(hrdetect_output[reorder,col_needed])
  
  par(mar=c(8,6,1,4))
  barplot(height = matrix_to_plot,
          beside = TRUE,
          #names.arg = row.names(hrdetect_res),
          names.arg = rownames(hrdetect_output)[reorder],
          las = 3,
          legend.text = c("deletion with MH",
                          "Substitution Sig. 3",
                          "Rearrangement Sig. 3",
                          "Rearrangement Sig. 5",
                          "HRD-LOH score",
                          "Substitution Sig. 8"),
          args.legend = list(x ='bottom', bty='n', inset=c(0,-0.9),horiz=TRUE,cex=1.5),
          ylim = c(min(matrix_to_plot)*1.1,max(matrix_to_plot)*1.1),
          col=c("blue","red","black","green","orange","yellow"),
          ylab = "BRCA1/BRCA2\n def. contribution",
          cex.axis = 1.2,
          cex.lab = 1.5,
          cex.names = 1)
  
  dev.off()
}

#----------------------

#' plot_HRDetect_overall
#' 
#' Overall plot of scores obtained from HRDetect.
#' 
#' @param file_name name of the output file (jpg)
#' @param hrdetect_output output of HRDetect, containing contributions of each feature
#' @export
plot_HRDetect_overall <- function(file_name,hrdetect_output){
  
  hrdetect_output <- as.data.frame(hrdetect_output)
  
  #barplot
  plot_colours <- rep("grey",nrow(hrdetect_output))
  
  reorder <- order(hrdetect_output$Probability,decreasing = TRUE)
  plot_colours <- plot_colours[reorder]
  jpeg(filename = file_name,
       width = 1500,
       height = 900,
       res = 200)
  par(xpd=FALSE)
  bp <- barplot(height = hrdetect_output$Probability[reorder],
                names.arg = "",
                #las = 2,
                main = paste0("HRDetect probability score"),
                ylim = c(0,1),
                border = 0,
                space = 0,
                col = plot_colours,ylab = "score")
  abline(h=0.7,col="red")
  par(xpd=TRUE)
  start <- 0
  for(i in 1:length(plot_colours)){
    rect(start, -0.1, start+1, -0.02,col = plot_colours[start+1],lwd = 0)
    start <- start + 1
  }
  # legend("bottom",c("BRCA1","BRCA2","BRCA1 Meth","BRCA2.note","BRCA1/MUTYH","PALB2"),inset=c(0,-0.2), horiz = TRUE, col=c("red","blue","green","purple","orange","yellow"),bty = "n", cex = 0.7,lwd=5)
  par(xpd=FALSE)
  dev.off()
}

#----------------------

#' plot_HRDetect_BootstrapScores
#' 
#' Overall plot of scores obtained from HRDetect with bootstrap.
#' 
#' @param outdir output directory for the plot
#' @param hrdetect_res output of HRDetect_pipeline, ran with bootstrap_scores=TRUE
#' @param main plot title. If not specified: "Distribution of HRDetect scores"
#' @param mar plot margin. If not specified: c(9,4,3,2)
#' @param pwidth plot width in pixels. If not specified: max(1000,400+120*number_of_samples)
#' @param pheight plot height in pixels. If not specified: 1000
#' @param pres plot resolution. If not specified: 200
#' @export
plot_HRDetect_BootstrapScores <- function(outdir,
                                          hrdetect_res,
                                          main="Distribution of HRDetect scores",
                                          mar=c(9,4,3,2),
                                          pwidth=NULL,
                                          pheight=1000,
                                          pres=200){
  
  #boxplot
  bootstrap_samples <- rownames(hrdetect_res$q_5_50_95)
  samples_labels <- bootstrap_samples
  
  if(is.null(pwidth)) pwidth <- max(1000,400+120*length(bootstrap_samples))
  
  o <- order(hrdetect_res$hrdetect_output[bootstrap_samples,"Probability"],decreasing = TRUE)
  jpeg(filename = paste0(outdir,"/HRDetect_bootstrap.jpg"),width = pwidth,height = pheight,res = pres)
  par(mar=mar)
  boxplot(hrdetect_res$hrdetect_bootstrap_table, xaxt="n", yaxt="n",ylim = c(0,1),border="white")
  grid()
  lines(x=c(0,ncol(hrdetect_res$hrdetect_bootstrap_table)+1),y=c(0.7,0.7),lty=2,col="red")
  boxplot(hrdetect_res$hrdetect_bootstrap_table[,o],ylim = c(0,1),main="",
          ylab="HRDetect score",las=3,names=samples_labels[o],add = TRUE)
  thisvalues <- hrdetect_res$hrdetect_bootstrap_table
  # take the x-axis indices and add a jitter, proportional to the N in each level
  for(i in 1:ncol(thisvalues)){
    myjitter<-jitter(rep(i, length(thisvalues[,o[i]])), amount=0.2)
    points(myjitter, thisvalues[,o[i]], pch=20, col=rgb(0,0,1,.1))
    rect(xleft = i+0.15,ybottom = hrdetect_res$q_5_50_95[o[i],"5%"],xright = i+0.22,ytop = hrdetect_res$q_5_50_95[o[i],"95%"],col = "grey")
  }
  points(hrdetect_res$hrdetect_output[bootstrap_samples,"Probability"][o],pch=23,col="black",bg="white")
  
  title(main = main,line = 1.8)
  legend(x="topleft",legend = c("Davies et al. 2017"),pch = 23,col="black",bg="green",cex = 0.8,bty = "n",inset = c(0,-0.145),xpd = TRUE)
  legend(x="topright",legend = c("5-95% CI"),col="black",fill="grey",cex = 0.8,bty = "n",inset = c(0,-0.145),xpd = TRUE)
  legend(x="topleft",legend = c("classification threshold, Davies et al. 2017"),lty = 2,col="red",cex = 0.8,bty = "n",inset = c(0,-0.095),xpd = TRUE)
  legend(x="topright",legend = c(paste0("n=",nrow(hrdetect_res$hrdetect_bootstrap_table))),lty = 2,col="white",cex = 0.8,bty = "n",inset = c(0,-0.095),xpd = TRUE)
  
  dev.off()
  
}

#' plot_HRDetect_Contributions
#' 
#' Function for plotting HRDetect results and contributions to HRDetect score.
#' 
#' @param file_name name of the output file (jpg)
#' @param par_title title of the plot (optional) 
#' @param hrdetect_output output of HRDetect, containing contributions of each feature
#' @export
plot_HRDetect_Contributions <- function(file_name,
                                                      hrdetect_output,
                                                      par_title = "HRDetect Score and Contributions"){
  
  col_needed <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
  
  #details plot
  hrdetect_output <- as.data.frame(hrdetect_output)
  reorder <- order(hrdetect_output$Probability,decreasing = TRUE)
  nsamples <- nrow(hrdetect_output)
  contributions <- hrdetect_output
  jpeg(filename = file_name,
       width = max(200+120*nsamples,2400),
       height = 1200,
       res = 170)
  par(mfrow=c(1,1),oma=c(0,0,1,0))
  mat <- matrix(c(1,2),ncol = 1)
  layout(mat, c(1), c(0.65,1))
  
  par(mar=c(1,6,1,1))
  barplot(height = hrdetect_output$Probability[reorder],
          ylim = c(0,1),ylab = "BRCA1/BRCA2\n def. score",
          names.arg = "",
          cex.axis = 1.2,
          cex.lab=1.2,
          cex.names = 1.5)
  title(main=paste0(par_title), cex.main=1.6)
  abline(a=0.7,b=0,col="red")
  text(x = 0.75*nrow(hrdetect_output),y = 0.77,labels = "classification threshold",col = "red",cex = 1.5)
  matrix_to_plot <- t(hrdetect_output[reorder,col_needed])
  
  par(mar=c(10,6,1,1))
  barplot(height = matrix_to_plot,
          beside = TRUE,
          #names.arg = row.names(hrdetect_res),
          names.arg = rownames(hrdetect_output)[reorder],
          las = 3,
          legend.text = c("deletion with MH",
                          "Substitution Sig. 3",
                          "Rearrangement Sig. 3",
                          "Rearrangement Sig. 5",
                          "HRD-LOH score",
                          "Substitution Sig. 8"),
          args.legend = list(x ='bottom', bty='n', inset=c(0,-1.0),horiz=TRUE,cex=0.9),
          # ylim = c(min(matrix_to_plot)*1.1,max(matrix_to_plot)*1.1),
          ylim = c(-3.5,4.5),
          col=c("blue","red","black","green","orange","yellow"),
          ylab = "BRCA1/BRCA2\n def. contribution",
          cex.axis = 1.2,
          cex.lab = 1.2,
          cex.names = 0.8)
  
  dev.off()
}

