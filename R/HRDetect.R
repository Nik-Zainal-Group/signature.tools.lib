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
#' will be returned as well. The input data_matrix can also be omitted, and the required
#' HRDetect features will be computed from the other input files. Signature fit with bootstrap
#' and bootstrap HRDetect scores can be requested using the bootstrapSignatureFit and
#' bootstrapHRDetectScores parameters.
#'
#' Single Nucleotide Variations. Columns in data_matrix relative to SNV are SNV3 and SNV8. Values
#' corresponding to number of SNV3 and SNV8 mutations in each sample can be provided in the data frame data_matrix.
#' Alternatively, SNV catalogues can be provided for the samples
#' (96-channels as rows and samples as columns) or can be constructed providing a list of either SNV VCF files or SNV TAB files.
#' The SNV catalogues will then be used to estimate signature exposures, using either RefSig or
#' COSMIC signatures. The signature version can be specified using the SNV_signature_version
#' parameter and a specific subset of signatures can be requested via the SNV_signature_names
#' parameter. If an organ is specified, the pipeline will attempt to use organ specific signatures.
#' If the signature version is RefSigv2 and an organ is specified, signature fit will be performed
#' using the FitMS algorithm (Degasperi et al. 2022, Science).
#'
#' Structural Variants (Rearrangements). Columns in data_matrix relative to SV are SV3 and SV5. Values
#' corresponding to number of SV3 and SV5 rearrangements in each sample can be provided in the data frame data_matrix.
#' Alternatively, SV catalogues can be provided for the samples
#' (32-channels as rows and samples as columns) or can be constructed providing a list of SV BEDPE files.
#' The SV catalogues will then be used to estimate signature exposures, using RefSig signatures.
#' The signature version can be specified using the SNV_signature_version
#' parameter and a specific subset of signatures can be requested via the SV_signature_names
#' parameter. If an organ is specified, the pipeline will attempt to use organ specific signatures.
#'
#' If signature fit for SNV or SV has already been performed using the Fit or FitMS functions,
#' the resulting objects can be passed directly using the subs_fit_obj and rearr_fit_obj
#' parameters. If the objects contain bootstrap fits, these will be used when the
#' bootstrap HRDetect score is requested. The HRDetect_pipeline function will attempt to extract
#' values for SNV3, SNV8, SV3, SV5, using the following signature names:
#' SNV3 = "SBS3", "Signature3", "RefSig3";
#' SNV8 = "SBS8", "Signature8", "RefSig8";
#' SV3 = "RS3","RefSigR3";
#' SV5 = "RS5", "RefSigR5", "RefSigR9".
#' If custom signature names have been used in the subs_fit_obj and rearr_fit_obj fits,
#' then the custom names can be provided using the parameters customNameSNV3, customNameSNV8,
#' customNameSV3 and customNameSV5.
#'
#' Deletions at Micro-homology (Indels). The column in data_matrix corresponding to the proportion of deletions at micro-homology is del.mh.prop.
#' The proportion of deletions at micro-homology for the samples can be calculated by the pipeline if the user provides Indels VCF files.
#'
#' HRD-LOH index (CNV). The column in data_matrix corresponding to the HRD-LOH index is hrd.
#' The HRD-LOH index for the samples can be calculated by the pipeline if the user provides copy numbers TAB files.
#'
#' The pipeline will produce some feedback in the form or info, warning, and error messages.
#' Please check the output to see whether everything worked as planned.
#'
#' @param data_matrix data frame containing a sample for each row and the six necessary features as columns. Columns should be labelled with the following names: del.mh.prop, SNV3, SV3, SV5, hrd, SNV8. Row names of the data frame should correspond to the sample names. If the values of the features need to be computed, set them to NA and provide additional data (e.g. catalogues, VCF/BEDPE/TAB files as specified in this documentation page).
#' @param genome.v genome version to use when constructing the SNV catalogue and classifying indels. Set it to either "hg19" or "hg38".
#' @param SNV_catalogues data frame containing 96-channel substitution catalogues. A sample for each column and the 96-channels as rows. Row names should have the correct channel names (see for example tests/testthat/test.snv.tab) and the column names should be the sample names so that each catalogue can be matched with the corresponding row in the data_matrix input.
#' @param SV_catalogues data frame containing 32-channel substitution catalogues. A sample for each column and the 32-channels as rows. Row names should have the correct channel names (see for example tests/testthat/test.cat) and the column names should be the sample names so that each catalogue can be matched with the corresponding row in the data_matrix input.
#' @param SNV_vcf_files list of file names corresponding to SNV VCF files to be used to construct 96-channel substitution catalogues. This should be a named vector, where the names indicate the sample name, so that each file can be matched to the corresponding row in the data_matrix input. The files should only contain SNV and should already be filtered according to the user preference, as all SNV in the file will be used and no filter will be applied.
#' @param SNV_tab_files list of file names corresponding to SNV TAB files to be used to construct 96-channel substitution catalogues. This should be a named vector, where the names indicate the sample name, so that each file can be matched to the corresponding row in the data_matrix input. The files should only contain SNV and should already be filtered according to the user preference, as all SNV in the file will be used and no filter will be applied. The files should contain a header in the first line with the following columns: chr, position, REF, ALT.
#' @param Indels_vcf_files list of file names corresponding to Indels VCF files to be used to classify Indels and compute the proportion of indels at micro-homology. This should be a named vector, where the names indicate the sample name, so that each file can be matched to the corresponding row in the data_matrix input. The files should only contain indels (no SNV) and should already be filtered according to the user preference, as all indels in the file will be used and no filter will be applied.
#' @param Indels_tab_files list of file names corresponding to Indels TAB files to be used to classify Indels and compute the proportion of indels at micro-homology. This should be a named vector, where the names indicate the sample name, so that each file can be matched to the corresponding row in the data_matrix input. The files should only contain indels (no SNV) and should already be filtered according to the user preference, as all indels in the file will be used and no filter will be applied. Each File contains indels from a single sample and the following minimal columns: chr, position, REF, ALT.
#' @param CNV_tab_files list of file names corresponding to CNV TAB files (similar to ASCAT format) to be used to compute the HRD-LOH index. This should be a named vector, where the names indicate the sample name, so that each file can be matched to the corresponding row in the data_matrix input. The files should contain a header in the first line with the following columns: 'seg_no', 'Chromosome', 'chromStart', 'chromEnd', 'total.copy.number.inNormal', 'minor.copy.number.inNormal', 'total.copy.number.inTumour', 'minor.copy.number.inTumour'
#' @param SV_bedpe_files list of file names corresponding to SV (Rearrangements) BEDPE files to be used to construct 32-channel rearrangement catalogues. This should be a named vector, where the names indicate the sample name, so that each file can be matched to the corresponding row in the data_matrix input. The files should contain a rearrangement for each row (two breakpoint positions should be on one row as determined by a pair of mates of paired-end sequencing) and should already be filtered according to the user preference, as all rearrangements in the file will be used and no filter will be applied. The files should contain a header in the first line with the following columns: "chrom1", "start1", "end1", "chrom2", "start2", "end2" and "sample" (sample name). In addition, either two columns indicating the strands of the mates, "strand1" (+ or -) and "strand2" (+ or -), or one column indicating the structural variant class, "svclass": translocation, inversion, deletion, tandem-duplication. The column "svclass" should correspond to (Sanger BRASS convention): inversion (strands +/- or -/+ and mates on the same chromosome), deletion (strands +/+ and mates on the same chromosome), tandem-duplication (strands -/- and mates on the same chromosome), translocation (mates are on different chromosomes)..
#' @param organ when using RefSigv1 or RefSigv2 as SNV_signature_version, organ-specific signatures will be used. Use one of the following organs: "Biliary", "Bladder", "Bone_SoftTissue", "Breast", "Cervix" (v1 only), "CNS", "Colorectal", "Esophagus", "Head_neck", "Kidney", "Liver", "Lung", "Lymphoid", "NET" (v2 only), "Oral_Oropharyngeal" (v2 only), "Ovary", "Pancreas", "Prostate", "Skin", "Stomach", "Uterus". If a certain organ is not available for either SNV or Sv signatures, or if the organ parameter is unspecified, then the pipeline will attempt to use the corresponding reference signatures instead, and SNV_signature_names and/or SV_signature_names can be used to specify a subset of signature names.
#' @param SNV_signature_version version of single base substitution signatures to use, either "COSMICv2", "COSMICv3.2", "RefSigv1" (Degasperi et al. 2020, Nature Cancer) or "RefSigv2" (Degasperi et al. 2022, Science)
#' @param SV_signature_version version of rearrangement signatures to use, only "RefSigv1" (Degasperi et al. 2020, Nature Cancer) currently available
#' @param SNV_signature_names when organ is not specified, you can use this to specify a list of SNV signature names to select from the set of signatures determined by the SNV_signature_version option
#' @param SV_signature_names when organ is not specified, you can use this to specify a list of SNV signature names to select from the set of signatures determined by the SV_signature_version option
#' @param subs_fit_obj Fit or FitMS result object. This parameter should be used when the user wants to customise the subs fit outside the HRDetect pipeline. If custom signatures were used, parameters customNameSNV3 and customNameSNV8 can be used to specify which custom signatures correspond to the HRDetect parameters SNV3 and SNV8.
#' @param rearr_fit_obj Fit or FitMS result object. This parameter should be used when the user wants to customise the rearrangements fit outside the HRDetect pipeline. If custom signatures were used, parameters customNameSV3 and customNameSV5 can be used to specify which custom signatures correspond to the HRDetect parameters SV3 and SV5.
#' @param customNameSNV3 custom signature name that will be considered as SNV3 input for HRDetect. Useful for when subs_fit_obj is provided and custom signatures are used.
#' @param customNameSNV8 custom signature name that will be considered as SNV8 input for HRDetect. Useful for when subs_fit_obj is provided and custom signatures are used.
#' @param customNameSV3 custom signature name that will be considered as SV3 input for HRDetect. Useful for when rearr_fit_obj is provided and custom signatures are used.
#' @param customNameSV5 custom signature name that will be considered as SV5 input for HRDetect. Useful for when rearr_fit_obj is provided and custom signatures are used.
#' @param optimisation_method can be KLD (KL divergence), NNLS (non-negative least squares) or SA (simulated annealing)
#' @param exposureFilterTypeFit use either fixedThreshold or giniScaledThreshold as exposure filter in signature fit. When using fixedThreshold, exposures will be removed based on a fixed percentage with respect to the total number of mutations (threshold_percentFit will be used). When using giniScaledThreshold each signature will used a different threshold calculated as (1-Gini(signature))*giniThresholdScalingFit
#' @param giniThresholdScalingFit scaling factor for when exposureFilterTypeFit="giniScaledThreshold", which is based on the Gini score of a signature
#' @param threshold_percentFit threshold in percentage of total mutations in a sample for when exposureFilterTypeFit="fixedThreshold". Only exposures larger than or equal to the threshold are considered, the others are set to zero
#' @param bootstrapSignatureFit set to TRUE to compute bootstrap signature fits, otherwise FALSE will compute a single fit. If a sample has a low number of mutations, then the bootstrap procedure can alter the catalogue a lot, in which case a single fit is advised
#' @param nbootFit number of bootstraps to use, more bootstraps more accurate results
#' @param threshold_p.valueFit p-value to determine whether an exposure is above the threshold_percent. In other words, this is the empirical probability that the exposure is lower than the threshold
#' @param bootstrapHRDetectScores perform HRDetect score with bootstrap. This requires mutations or catalogues for subs/rearr to compute the bootstrap fit, and indels mutations to bootstrap the indels classification. HRD-LOH can still be provided using the input data_matrix.
#' @param nparallel how many parallel threads to use.
#' @param randomSeed set an integer random seed
#' @return return a list that contains $data_matrix (updated input data_matrix with additional computed features), $hrdetect_output (data frame with HRDetect BRCAness Probability and contribution of the features), $SNV_catalogues (input SNV_catalogues updated with additional computed substitution catalogues if any), $SV_catalogues (input SV_catalogues updated with additional computed rearrangement catalogues if any)
#' @export
#' @references A. Degasperi, T. D. Amarante, J. Czarnecki, S. Shooter, X. Zou, D. Glodzik, ... H. Davies, S. Nik-Zainal. A practical framework and online tool for mutational signature analyses show intertissue variation and driver dependencies, Nature Cancer, https://doi.org/10.1038/s43018-020-0027-5, 2020.
#' @references A. Degasperi, X. Zou, T. D. Amarante, ..., H. Davies, Genomics England Research Consortium, S. Nik-Zainal. Substitution mutational signatures in whole-genome-sequenced cancers in the UK population. Science, 2022.
#' @references Davies, H., Glodzik, D., Morganella, S., Yates, L. R., Staaf, J., Zou, X., ... Nik-Zainal, S. (2017). HRDetect is a predictor of BRCA1 and BRCA2 deficiency based on mutational signatures. Nature Medicine, 23(4), 517–525. https://doi.org/10.1038/nm.4292
#' @references Nik-Zainal, S., Davies, H., Staaf, J., Ramakrishna, M., Glodzik, D., Zou, X., ... Stratton, M. R. (2016). Landscape of somatic mutations in 560 breast cancer whole-genome sequences. Nature, 534(7605), 1–20. https://doi.org/10.1038/nature17676
#' @references Huang, X., Wojtowicz, D., & Przytycka, T. M. (2017). Detecting Presence Of Mutational Signatures In Cancer With Confidence. bioRxiv, (October). https://doi.org/10.1101/132597
#'
HRDetect_pipeline <- function(data_matrix=NULL,
                              genome.v="hg19",
                              SNV_catalogues=NULL,
                              SV_catalogues=NULL,
                              SNV_vcf_files=NULL,
                              SNV_tab_files=NULL,
                              Indels_vcf_files=NULL,
                              Indels_tab_files=NULL,
                              CNV_tab_files=NULL,
                              SV_bedpe_files=NULL,
                              organ=NULL,
                              SNV_signature_version="RefSigv2",
                              SV_signature_version="RefSigv1",
                              SNV_signature_names=NULL,
                              SV_signature_names=NULL,
                              subs_fit_obj=NULL,
                              rearr_fit_obj=NULL,
                              customNameSNV3=NULL,
                              customNameSNV8=NULL,
                              customNameSV3=NULL,
                              customNameSV5=NULL,
                              optimisation_method = "KLD",
                              exposureFilterTypeFit = "fixedThreshold",
                              giniThresholdScalingFit = 10,
                              threshold_percentFit = 5,
                              bootstrapSignatureFit = TRUE,
                              nbootFit=100,
                              threshold_p.valueFit = 0.05,
                              bootstrapHRDetectScores=FALSE,
                              nparallel=1,
                              randomSeed=NULL){
  #if multiple parallel cores are used, set it here
  doParallel::registerDoParallel(nparallel)

  if(!is.null(randomSeed)){
    set.seed(randomSeed)
  }

  message("[info HRDetect_pipeline] HRDetect pipeline starting!")

  # check for signature fit files provided
  custom_subsFit <- NULL
  custom_rearrFit <- NULL
  if(!is.null(subs_fit_obj)){
    if(!is.null(subs_fit_obj$exposures) & !is.null(subs_fit_obj$fitAlgorithm)){
      message("[info HRDetect_pipeline] using subs_fit_obj ",subs_fit_obj$fitAlgorithm," object provided.")
      custom_subsFit <- subs_fit_obj
    }else{
      message("[error HRDetect_pipeline] subs_fit_obj object provided does not appear to be a valid Fit or FitMS object.")
      return(NULL)
    }

  }
  if(!is.null(rearr_fit_obj)){
    if(!is.null(rearr_fit_obj$exposures) & !is.null(rearr_fit_obj$fitAlgorithm)){
      message("[info HRDetect_pipeline] using rearr_fit_obj ",rearr_fit_obj$fitAlgorithm," object provided.")
      custom_rearrFit <- rearr_fit_obj
    }else{
      message("[error HRDetect_pipeline] rearr_fit_obj object provided does not appear to be a valid Fit or FitMS object.")
      return(NULL)
    }
  }

  #check that the matrix has correct features (columns)
  col_needed <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
  # need to check if there are samples that require fit
  SNV_samples_with_data <- unique(unlist(sapply(list(SNV_vcf_files,
                                                     SNV_tab_files,
                                                     SNV_catalogues), function(x) names(x))))
  SV_samples_with_data <- unique(unlist(sapply(list(SV_bedpe_files,
                                                    SV_catalogues), function(x) names(x))))
  samples_with_data <- unique(unlist(sapply(list(SNV_vcf_files,
                                                    SNV_tab_files,
                                                    SNV_catalogues,
                                                    Indels_vcf_files,
                                                    Indels_tab_files,
                                                    CNV_tab_files,
                                                    SV_bedpe_files,
                                                    SV_catalogues), function(x) names(x))))
  if(!is.null(custom_subsFit)) samples_with_data <- union(samples_with_data,rownames(custom_subsFit$exposures))
  if(!is.null(custom_rearrFit)) samples_with_data <- union(samples_with_data,rownames(custom_rearrFit$exposures))

  if (!is.null(data_matrix)){
    if (!length(intersect(col_needed,colnames(data_matrix)))==length(col_needed)){
      message("[error HRDetect_pipeline] incorrect data_matrix columns specified, you need the following columns: \"del.mh.prop\", \"SNV3\", \"SV3\", \"SV5\", \"hrd\", \"SNV8\"")
      return(NULL)
    }
    # check that all samples with data are in the rownames of the data_matrix
    sample_names_missing <- setdiff(samples_with_data,rownames(data_matrix))
    if(length(sample_names_missing)>0){
      message("[warning HRDetect_pipeline] some sample names associated with mutation/catalogue files are missing in the data_matrix provided. These will be added to the data_matrix as additional rows.")
      data_matrix <- rbind(data_matrix,matrix(NA,ncol = 6,nrow = length(sample_names_missing),dimnames = list(sample_names_missing,col_needed)))
    }
  }else{
    #there is no data_matrix, just build it from files
    if(!is.null(samples_with_data)){
      data_matrix <- matrix(NA,ncol = 6,nrow = length(samples_with_data),dimnames = list(samples_with_data,col_needed))
    }else{
      message("[error HRDetect_pipeline] no input data provided.")
      return(NULL)
    }
  }
  samples_list <- rownames(data_matrix)

  #check that the signature version options are valid

  if(!(SNV_signature_version %in% c("COSMICv2","COSMICv3.2","RefSigv1","RefSigv2"))){
    message(paste0("[error HRDetect_pipeline] invalid SNV_signature_version ",SNV_signature_version,"."))
    return(NULL)
  }
  if(!(SV_signature_version %in% c("RefSigv1"))){
    message(paste0("[error HRDetect_pipeline] invalid SV_signature_version ",SV_signature_version,"."))
    return(NULL)
  }

  #initialise bootstrap fit variables to NULL
  fitRes_subs <- NULL
  fitRes_rearr <- NULL
  #initialise final exposures output table
  exposures_subs <- NULL
  exposures_rearr <- NULL
  #initialise indels classification table
  indels_classification_table <- NULL
  #initialise hrdetect bootstrap tables
  hrdetect_bootstrap_table <- NULL
  q_5_50_95 <- NULL
  # annotated mutations
  annotated_mutations_subs <- NULL
  annotated_mutations_rearr <- NULL

  # first of all, we need to check that signature fit is necessary
  SNV_cols <- c("SNV3","SNV8")
  SV_cols <- c("SV3","SV5")
  SNV3names <- c("SBS3","Signature3","RefSig3")
  SNV8names <- c("SBS8","Signature8","RefSig8")
  SV3names <- c("RS3","RefSigR3")
  SV5names <- c("RS5","RefSigR5","RefSigR9")
  # add custom names for the above HRDetect input signatures
  if(!is.null(customNameSNV3)) SNV3names <- union(SNV3names,customNameSNV3)
  if(!is.null(customNameSNV8)) SNV8names <- union(SNV8names,customNameSNV8)
  if(!is.null(customNameSV3)) SV3names <- union(SV3names,customNameSV3)
  if(!is.null(customNameSV5)) SV5names <- union(SV5names,customNameSV5)

  message("[info HRDetect_pipeline] Single Nucleotide Variations")

  # check if we have custom subs signature fit
  if(!is.null(custom_subsFit)){
    # check for overlapping sample names, custom fit has the precedence
    custom_subs_samples <- rownames(custom_subsFit$exposures)
    subs_conflict_samples <- intersect(SNV_samples_with_data,custom_subs_samples)
    if(length(subs_conflict_samples)>0){
      message("[warning HRDetect_pipeline] There are sample names conflicts for substitutions. Some samples provided via subs_fit_file have also been provided via either catalogues or mutation files. Results in subs_fit_file will be used for samples ",paste(subs_conflict_samples,collapse = ","),".")
      SNV_samples_with_data <- setdiff(SNV_samples_with_data,subs_conflict_samples)
    }
  }

  # select appropriate SNV signatures to use and run signatureFit_pipeline
  if(!is.null(SNV_samples_with_data)){

    SNV_fit_method <- "Fit"

    if(SNV_signature_version=="RefSigv2"){
      if(!is.null(organ)){
        # check if the organ is available in FitMS
        available_organs <- rownames(sigsForFittingSBSv2.03)
        if(organ %in% available_organs){
          SNV_fit_method <- "FitMS"
        }else{
          message("[warning HRDetect_pipeline] SNV organ-specific signatures not available for ",organ,", version RefSigv2, trying to fix this by switching to use the reference signatures instead. Please use SNV_signature_names to specify a subset of reference signatures to use.")
          organ <- NULL
        }
      }
    }else if(SNV_signature_version=="RefSigv1"){
      if(!is.null(organ)){
        # check if the organ is available
        organsigs <- getOrganSignatures(organ = organ,version = "1",typemut = "subs",verbose = FALSE)
        if(ncol(organsigs)==0){
          message("[warning HRDetect_pipeline] SNV organ-specific signatures not available for ",organ,", version RefSigv1, trying to fix this by switching to use the reference signatures instead. Please use SNV_signature_names to specify a subset of reference signatures to use.")
          organ <- NULL
        }
      }
    }
    # if SNV_signature_version is COSMICv2 or COSMICv3.2, there is currently nothing to do
    # the function signatureFit_pipeline will take care of it, and if organ is not NULL
    # we can add some checks when the organ option is implemented for COSMIC signatures

    message("[info HRDetect_pipeline] Running SNV signature extraction using the ",SNV_fit_method," method for samples ",paste(SNV_samples_with_data,collapse = ", "),".")
    # Now I can run the signature fit using signatureFit_pipeline
    fitPipeline_SNV <- signatureFit_pipeline(catalogues = SNV_catalogues,
                                             genome.v = genome.v,
                                             organ = organ,
                                             SNV_vcf_files = SNV_vcf_files,
                                             SNV_tab_files = SNV_tab_files,
                                             signature_version = SNV_signature_version,
                                             signature_names = SNV_signature_names,
                                             fit_method = SNV_fit_method,
                                             optimisation_method = optimisation_method,
                                             useBootstrap = bootstrapSignatureFit,
                                             exposureFilterType = exposureFilterTypeFit,
                                             threshold_percent = threshold_percentFit,
                                             giniThresholdScaling = giniThresholdScalingFit,
                                             nboot = nbootFit,
                                             nparallel = nparallel,
                                             randomSeed = randomSeed)

    if(is.null(fitPipeline_SNV$fitResults)){
      message("[error HRDetect_pipeline] SNV signatures fit failed.")
      return(NULL)
    }else{
      message("[info HRDetect_pipeline] SNV signatures fit completed.")
      # save the results
      SNV_catalogues <- fitPipeline_SNV$catalogues
      fitRes_subs <- fitPipeline_SNV$fitResults
      exposures_subs <- t(fitRes_subs$exposures)
      exposures_subs <- exposures_subs[1:(nrow(exposures_subs)-1),,drop=F]
      annotated_mutations_subs <- fitPipeline_SNV$annotated_mutations
    }

  }else{
    message("[info HRDetect_pipeline] no SNV catalogues or SNV mutations provided, no signature fit performed.")
  }

  # if we have new exposures we can add them to the data_matrix
  if(!is.null(exposures_subs)){

    # set default value
    data_matrix[colnames(exposures_subs),"SNV3"] <- 0
    data_matrix[colnames(exposures_subs),"SNV8"] <- 0

    # if SNV RefSigv1 or RefSigv2 are used AND organ is not null, then we need to convert to reference signatures
    if((SNV_signature_version=="RefSigv1" | SNV_signature_version=="RefSigv2") & !is.null(organ)){
      # if RefSigv2 was used, there may be some SBS sigs in the names
      sbssigs <- rownames(exposures_subs)[grepl(rownames(exposures_subs),pattern = "^SBS")]
      if(length(sbssigs)>0){
        exposures_subs_toconvert <- exposures_subs[setdiff(rownames(exposures_subs),sbssigs),,drop=F]
        exposures_subs_sbs <- exposures_subs[sbssigs,,drop=F]
        exposures_subs_converted <- convertExposuresFromOrganToRefSigs(exposures_subs_toconvert,typemut = "subs")
        exposures_subs_converted <- exposures_subs_converted[apply(exposures_subs_converted, 1,sum)>0,,drop=FALSE]
        exposures_subs <- rbind(exposures_subs_converted,exposures_subs_sbs)
      }else{
        exposures_subs <- convertExposuresFromOrganToRefSigs(exposures_subs,typemut = "subs")
        exposures_subs <- exposures_subs[apply(exposures_subs, 1,sum)>0,,drop=FALSE]
      }
    }

    # update the data_matrix
    data_matrix <- updateDataMatrix(data_matrix = data_matrix,
                                    exposures = exposures_subs,
                                    data_matrixFeature = "SNV3",
                                    SigNames = SNV3names)
    data_matrix <- updateDataMatrix(data_matrix = data_matrix,
                                    exposures = exposures_subs,
                                    data_matrixFeature = "SNV8",
                                    SigNames = SNV8names)
  }

  # check for custom SNV signature fit
  if(!is.null(custom_subsFit)){

    custom_exposures_subs <- t(custom_subsFit$exposures)
    custom_exposures_subs <- custom_exposures_subs[1:(nrow(custom_exposures_subs)-1),,drop=F]

    # update data matrix
    data_matrix <- updateDataMatrix(data_matrix = data_matrix,
                                    exposures = custom_exposures_subs,
                                    data_matrixFeature = "SNV3",
                                    SigNames = SNV3names)
    data_matrix <- updateDataMatrix(data_matrix = data_matrix,
                                    exposures = custom_exposures_subs,
                                    data_matrixFeature = "SNV8",
                                    SigNames = SNV8names)
  }


  message("[info HRDetect_pipeline] Structural Variants (Rearrangements)")

  # check if we have custom rearr signature fit
  if(!is.null(custom_rearrFit)){
    # check for overlapping sample names, custom fit has the precedence
    custom_rearr_samples <- rownames(custom_rearrFit$exposures)
    rearr_conflict_samples <- intersect(SV_samples_with_data,custom_rearr_samples)
    if(length(rearr_conflict_samples)>0){
      message("[warning HRDetect_pipeline] There are sample names conflicts for rearrangements. Some samples provided via rearr_fit_file have also been provided via either catalogues or mutation files. Results in rearr_fit_file will be used for samples ",paste(rearr_conflict_samples,collapse = ","),".")
      SV_samples_with_data <- setdiff(SV_samples_with_data,rearr_conflict_samples)
    }
  }

  # select appropriate SV signatures to use
  if(!is.null(SV_samples_with_data)){

    SV_fit_method <- "Fit"

    if(SV_signature_version=="RefSigv1"){
      if(!is.null(organ)){
        # check if the organ is available
        organsigs <- getOrganSignatures(organ = organ,version = "1",typemut = "rearr",verbose = FALSE)
        if(ncol(organsigs)==0){
          message("[warning HRDetect_pipeline] SV organ-specific signatures not available for ",organ,", version RefSigv1, trying to fix this by switching to use the reference signatures instead. Please use SV_signature_names to specify a subset of reference signatures to use.")
          organ <- NULL
        }
      }
    }

    message("[info HRDetect_pipeline] Running SV signature extraction using the ",SV_fit_method," method for samples ",paste(SV_samples_with_data,collapse = ", "),".")
    # Now I can run the signature fit using signatureFit_pipeline
    fitPipeline_SV <- signatureFit_pipeline(catalogues = SV_catalogues,
                                            genome.v = genome.v,
                                            organ = organ,
                                            SV_bedpe_files = SV_bedpe_files,
                                            signature_version = SV_signature_version,
                                            signature_names = SV_signature_names,
                                            fit_method = SV_fit_method,
                                            optimisation_method = optimisation_method,
                                            useBootstrap = bootstrapSignatureFit,
                                            exposureFilterType = exposureFilterTypeFit,
                                            threshold_percent = threshold_percentFit,
                                            giniThresholdScaling = giniThresholdScalingFit,
                                            nboot = nbootFit,
                                            nparallel = nparallel,
                                            randomSeed = randomSeed)

    if(is.null(fitPipeline_SV$fitResults)){
      message("[error HRDetect_pipeline] SV signatures fit failed.")
      return(NULL)
    }else{
      message("[info HRDetect_pipeline] SV signatures fit completed.")
      # save the results
      SV_catalogues <- fitPipeline_SV$catalogues
      fitRes_rearr <- fitPipeline_SV$fitResults
      exposures_rearr <- t(fitRes_rearr$exposures)
      exposures_rearr <- exposures_rearr[1:(nrow(exposures_rearr)-1),,drop=F]
      annotated_mutations_rearr <- fitPipeline_SV$annotated_mutations
    }

  }else{
    message("[info HRDetect_pipeline] no SV catalogues or SV mutations provided, no signature fit performed.")
  }

  # if we have new exposures we can add them to the data_matrix
  if(!is.null(exposures_rearr)){

    # set default
    data_matrix[colnames(exposures_rearr),"SV3"] <- 0
    data_matrix[colnames(exposures_rearr),"SV5"] <- 0

    # if SV RefSigv1 are used AND organ is not null, then we need to convert to reference signatures
    if((SV_signature_version=="RefSigv1") & !is.null(organ)){
      #convert to reference signatures
      exposures_rearr <- convertExposuresFromOrganToRefSigs(exposures_rearr,typemut = "rearr")
      exposures_rearr <- exposures_rearr[apply(exposures_rearr, 1,sum)>0,,drop=FALSE]
    }

    # update the data_matrix
    data_matrix <- updateDataMatrix(data_matrix = data_matrix,
                                    exposures = exposures_rearr,
                                    data_matrixFeature = "SV3",
                                    SigNames = SV3names)
    data_matrix <- updateDataMatrix(data_matrix = data_matrix,
                                    exposures = exposures_rearr,
                                    data_matrixFeature = "SV5",
                                    SigNames = SV5names)
  }

  # check for custom Rearr signature fit
  if(!is.null(custom_rearrFit)){

    custom_exposures_rearr <- t(custom_rearrFit$exposures)
    custom_exposures_rearr <- custom_exposures_rearr[1:(nrow(custom_exposures_rearr)-1),,drop=F]

    # update data matrix
    data_matrix <- updateDataMatrix(data_matrix = data_matrix,
                                    exposures = custom_exposures_rearr,
                                    data_matrixFeature = "SV3",
                                    SigNames = SV3names)
    data_matrix <- updateDataMatrix(data_matrix = data_matrix,
                                    exposures = custom_exposures_rearr,
                                    data_matrixFeature = "SV5",
                                    SigNames = SV5names)
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

    if (!is.null(Indels_tab_files)){
      incomplete_samples_with_tabIndels <- intersect(incomplete_samples,names(Indels_tab_files))
    }else{
      #there is no Indels tab files given, so no incomplete sample has a tab file
      incomplete_samples_with_tabIndels <- character(0)
    }

    #also if a sample has both vcf and tab file, use the vcf file
    incomplete_samples_with_tabIndels <- setdiff(incomplete_samples_with_tabIndels,incomplete_samples_with_vcfIndels)

    #compute del.mh.prop for the samples with vcf Indels file that do not have del.mh.prop yet
    if(length(incomplete_samples_with_vcfIndels)>0){
      message("[info HRDetect_pipeline] Proportion of Indels with MH will be computed for the following VCF samples: ",paste(incomplete_samples_with_vcfIndels,collapse = " "))

      mh_list <- foreach::foreach(sample=incomplete_samples_with_vcfIndels) %dopar% {
        res <- vcfToIndelsClassification(Indels_vcf_files[sample],sample,genome.v = genome.v)
        res$count_proportion
      }
      #combine in one table and add to data_matrix
      indels_classification_table <- do.call(rbind,mh_list)
      rownames(indels_classification_table) <- indels_classification_table[,"sample"]
      data_matrix[rownames(indels_classification_table),MH_cols] <- indels_classification_table[,"del.mh.prop"]
    }

    #compute del.mh.prop for the samples with tab Indels file that do not have del.mh.prop yet
    if(length(incomplete_samples_with_tabIndels)>0){
      message("[info HRDetect_pipeline] Proportion of Indels with MH will be computed for the following TAB samples: ",paste(incomplete_samples_with_tabIndels,collapse = " "))

      mh_list <- foreach::foreach(sample=incomplete_samples_with_tabIndels) %dopar% {
        indels_table <- read.table(Indels_tab_files[sample],sep = "\t",header = TRUE,stringsAsFactors = FALSE)
        res <- tabToIndelsClassification(indels_table,sample,genome.v = genome.v)
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
    message("[info HRDetect_pipeline] Some samples do not have data necessary to compute the HRDetect score (check output $data_matrix). Will NOT compute HRDetect score for the following samples: ",paste(incomplete_cases,collapse = " "))
    hrdetect_input <- data_matrix[complete.cases(data_matrix),,drop=FALSE]
  }else{
    hrdetect_input <- data_matrix
  }

  if(nrow(hrdetect_input)>0){
    message("[info HRDetect_pipeline] Computing HRDetect score and feature contributions for the following samples: ",paste(rownames(hrdetect_input),collapse = " "))
    hrdetect_output <- applyHRDetectDavies2017(hrdetect_input, attachContributions = TRUE)

    #attempt to run HRDetect with bootstrap if requested
    if(bootstrapHRDetectScores){
      message("[info HRDetect_pipeline] HRDetect boostrap scores requested!")

      #check whether we have all we need and for which samples
      bootstrap_samples <- c()

      #bootstrap exposures are required
      bootstrap_subs_samples <- c()
      if(!is.null(fitRes_subs$bootstrap_exposures_samples)) bootstrap_subs_samples <- names(fitRes_subs$bootstrap_exposures_samples)
      if(!is.null(custom_subsFit$bootstrap_exposures_samples)) bootstrap_subs_samples <- union(bootstrap_subs_samples,names(custom_subsFit$bootstrap_exposures_samples))
      bootstrap_rearr_samples <- c()
      if(!is.null(fitRes_rearr$bootstrap_exposures_samples)) bootstrap_rearr_samples <- names(fitRes_rearr$bootstrap_exposures_samples)
      if(!is.null(custom_rearrFit$bootstrap_exposures_samples)) bootstrap_rearr_samples <- union(bootstrap_rearr_samples,names(custom_rearrFit$bootstrap_exposures_samples))

      bootstrap_samples <- intersect(bootstrap_subs_samples,bootstrap_rearr_samples)

      #indels classification tables are required
      if(!is.null(indels_classification_table)) bootstrap_samples <- intersect(bootstrap_samples,rownames(indels_classification_table))

      #If samples are in the hrdetect_input table then they also have the HRD-LOH score
      bootstrap_samples <- intersect(bootstrap_samples,rownames(hrdetect_input))

      if(length(bootstrap_samples)>0){
        message("[info HRDetect_pipeline] Computing HRDetect boostrap for the following samples: ",paste(bootstrap_samples,collapse = ", "))

        #run hrdetect bootstrap scores for the bootstrap_samples
        nboots_hr <- 1000

        # set up some variables
        subs_boots_list <- list()
        subs_catalogues_list <- list()
        if(!is.null(fitRes_subs)){
          subs_boots_list[["fit_subs"]] <- fitRes_subs$bootstrap_exposures
          subs_catalogues_list[["fit_subs"]] <- fitRes_subs$catalogues[,intersect(bootstrap_samples,colnames(fitRes_subs$catalogues)),drop=F]
        }
        if(!is.null(custom_subsFit)){
          subs_boots_list[["custom_fit_subs"]] <- custom_subsFit$bootstrap_exposures
          subs_catalogues_list[["custom_fit_subs"]] <- custom_subsFit$catalogues[,intersect(bootstrap_samples,colnames(custom_subsFit$catalogues)),drop=F]
        }
        rearr_boots_list <- list()
        rearr_catalogues_list <- list()
        if(!is.null(fitRes_rearr)){
          rearr_boots_list[["fit_rearr"]] <- fitRes_rearr$bootstrap_exposures
          rearr_catalogues_list[["fit_rearr"]] <- fitRes_rearr$catalogues[,intersect(bootstrap_samples,colnames(fitRes_rearr$catalogue)),drop=F]
        }
        if(!is.null(custom_rearrFit)){
          rearr_boots_list[["custom_fit_rearr"]] <- custom_rearrFit$bootstrap_exposures
          rearr_catalogues_list[["custom_fit_rearr"]] <- custom_rearrFit$catalogues[,intersect(bootstrap_samples,colnames(custom_rearrFit$catalogue)),drop=F]
        }


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

          # there may be one or two bootstrap objects
          for(bn in names(subs_boots_list)){
            subs_boots <- subs_boots_list[[bn]]
            subs_catalogues <- subs_catalogues_list[[bn]]
            current_nbootFit <- length(subs_boots)

            s1 <- sample(current_nbootFit,size = 1)
            current_subs <- subs_boots[[s1]][,colnames(subs_catalogues),drop=F]

            # sparsity correction
            sums_exp <- apply(subs_catalogues, 2, sum)
            denominator <- matrix(sums_exp,nrow = nrow(current_subs),ncol = ncol(current_subs),byrow = TRUE)
            exposuresProp <- (current_subs/denominator*100)
            # case of empty catalogues
            exposuresProp[,sums_exp==0] <- 0
            current_subs[exposuresProp<threshold_percentFit] <- 0

            # set default value
            tmp_data_matrix[colnames(current_subs),"SNV3"] <- 0
            tmp_data_matrix[colnames(current_subs),"SNV8"] <- 0

            # if SNV RefSigv1 or RefSigv2 are used AND organ is not null, then we need to convert to reference signatures
            if((SNV_signature_version=="RefSigv1" | SNV_signature_version=="RefSigv2") & !is.null(organ)){
              # if RefSigv2 was used, there may be some SBS sigs in the names
              sbssigs <- rownames(current_subs)[grepl(rownames(current_subs),pattern = "^SBS")]
              if(length(sbssigs)>0){
                exposures_subs_toconvert <- current_subs[setdiff(rownames(current_subs),sbssigs),,drop=F]
                exposures_subs_sbs <- current_subs[sbssigs,,drop=F]
                exposures_subs_converted <- convertExposuresFromOrganToRefSigs(exposures_subs_toconvert,typemut = "subs")
                exposures_subs_converted <- exposures_subs_converted[apply(exposures_subs_converted, 1,sum)>0,,drop=FALSE]
                current_subs <- rbind(exposures_subs_converted,exposures_subs_sbs)
              }else{
                current_subs <- convertExposuresFromOrganToRefSigs(current_subs,typemut = "subs")
                current_subs <- current_subs[apply(current_subs, 1,sum)>0,,drop=FALSE]
              }
            }

            # update data matrix
            tmp_data_matrix <- updateDataMatrix(data_matrix = tmp_data_matrix,
                                                exposures = current_subs,
                                                data_matrixFeature = "SNV3",
                                                SigNames = SNV3names)
            tmp_data_matrix <- updateDataMatrix(data_matrix = tmp_data_matrix,
                                                exposures = current_subs,
                                                data_matrixFeature = "SNV8",
                                                SigNames = SNV8names)

          }

          #4) sample rearr exposures

          # there may be one or two bootstrap objects
          for(bn in names(rearr_boots_list)){
            rearr_boots <- rearr_boots_list[[bn]]
            rearr_catalogues <- rearr_catalogues_list[[bn]]
            current_nbootFit <- length(rearr_boots)

            s1 <- sample(current_nbootFit,size = 1)
            current_rearr <- rearr_boots[[s1]][,colnames(rearr_catalogues),drop=F]
            current_rearr[is.na(current_rearr)] <- 0
            # sparsity correction
            sums_exp <- apply(rearr_catalogues, 2, sum)
            denominator <- matrix(sums_exp,nrow = nrow(current_rearr),ncol = ncol(current_rearr),byrow = TRUE)
            exposuresProp <- (current_rearr/denominator*100)
            # # case of empty catalogues
            exposuresProp[,sums_exp==0] <- 0
            current_rearr[exposuresProp<threshold_percentFit] <- 0

            # set default
            tmp_data_matrix[colnames(current_rearr),"SV3"] <- 0
            tmp_data_matrix[colnames(current_rearr),"SV5"] <- 0

            # if SV RefSigv1 are used AND organ is not null, then we need to convert to reference signatures
            if((SV_signature_version=="RefSigv1") & !is.null(organ)){
              #convert to reference signatures
              current_rearr <- convertExposuresFromOrganToRefSigs(current_rearr,typemut = "rearr")
              current_rearr <- current_rearr[apply(current_rearr, 1,sum)>0,,drop=FALSE]
            }

            # update data matrix
            tmp_data_matrix <- updateDataMatrix(data_matrix = tmp_data_matrix,
                                                exposures = current_rearr,
                                                data_matrixFeature = "SV3",
                                                SigNames = SV3names)
            tmp_data_matrix <- updateDataMatrix(data_matrix = tmp_data_matrix,
                                                exposures = current_rearr,
                                                data_matrixFeature = "SV5",
                                                SigNames = SV5names)

          }
          # compute score
          hrdetect_bootstrap_table[[j]] <- t(applyHRDetectDavies2017(data_matrix = tmp_data_matrix,attachContributions = FALSE))
        }

        hrdetect_bootstrap_table <- do.call(rbind,hrdetect_bootstrap_table)
        rownames(hrdetect_bootstrap_table) <- 1:nrow(hrdetect_bootstrap_table)
        q_5_50_95 <- t(apply(hrdetect_bootstrap_table,2,function(x) quantile(x,c(0.05,0.5,0.95))))

        message("[info HRDetect_pipeline] HRDetect with bootstrap successful!")

      }else{
        message("[info HRDetect_pipeline] Not enough data to run HRDetect bootstrap scores.")
        if(!bootstrapSignatureFit) message("[info HRDetect_pipeline] bootstrapSignatureFit needs to be TRUE to compute HRDetect bootstrap scores.")
      }
    }

    message("[info HRDetect_pipeline] HRDetect pipeline completed!")
  }else{
    message("[info HRDetect_pipeline] Impossible to compute any HRDetect score, as no sample seems to have all the necessary data. Check output $data_matrix to see what has been computed with the data provided.")
    hrdetect_output <- NULL
  }


  #--- return results ---

  res <- list()
  res$parameters <- list()
  res$parameters$SNV_signature_version <- SNV_signature_version
  res$parameters$SV_signature_version <- SV_signature_version
  res$parameters$organ <- organ
  res$parameters$SNV_signature_names <- SNV_signature_names
  res$parameters$SV_signature_names <- SV_signature_names
  res$parameters$optimisation_method <- optimisation_method
  res$parameters$threshold_percentFit <- threshold_percentFit
  res$parameters$bootstrapSignatureFit <- bootstrapSignatureFit
  res$parameters$nbootFit <- nbootFit
  res$parameters$threshold_p.valueFit <- threshold_p.valueFit
  res$parameters$bootstrapHRDetectScores <- bootstrapHRDetectScores
  res$parameters$nparallel <- nparallel
  res$data_matrix <- data_matrix
  res$SV_catalogues <- SV_catalogues
  res$SNV_catalogues <- SNV_catalogues
  res$hrdetect_output <- hrdetect_output
  res$exposures_subs <- exposures_subs
  res$exposures_rearr <- exposures_rearr
  res$fitRes_subs <- fitRes_subs
  res$fitRes_rearr <- fitRes_rearr
  res$annotated_mutations_subs <- annotated_mutations_subs
  res$annotated_mutations_rearr <- annotated_mutations_rearr
  res$indels_classification_table <- indels_classification_table
  res$hrdetect_bootstrap_table <- hrdetect_bootstrap_table
  res$q_5_50_95 <- q_5_50_95
  return(res)

}

# function to update the data_matrix, to avoid code duplication
updateDataMatrix <- function(data_matrix,
                             exposures, #signatures are in the rows and samples as columns
                             data_matrixFeature, #SNV3, SNV8, SV5 or SV8
                             SigNames){ #signature names for the given feature
  # check for a given matrix column
  whichFeature <- rownames(exposures) %in% SigNames
  if(sum(whichFeature)>0){
    data_matrix[colnames(exposures),data_matrixFeature] <- apply(exposures[whichFeature,,drop=F],2,sum)
  }
  return(data_matrix)
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
#' Overall plot of scores obtained from HRDetect with bootstrap. Plot will be in pdf format.
#'
#' @param outdir output directory for the plot
#' @param hrdetect_res output of HRDetect_pipeline, ran with bootstrapHRDetectScores=TRUE
#' @param main plot title. If not specified: "Distribution of HRDetect scores"
#' @param mar plot margin. If not specified, this will be dynamically estimated based on sample names length
#' @param pwidth plot width in pixels. If not specified, this will be dynamically estimated based on the number of samples
#' @param pheight plot height in pixels. If not specified, this will be dynamically estimated based on sample names length
#' @param pointsize plot resolution. If not specified: 12
#' @export
plot_HRDetect_BootstrapScores <- function(outdir,
                                          hrdetect_res,
                                          main="Distribution of HRDetect scores",
                                          mar=NULL,
                                          pwidth=NULL,
                                          pheight=NULL,
                                          pointsize=12){

  #boxplot
  bootstrap_samples <- rownames(hrdetect_res$q_5_50_95)
  samples_labels <- bootstrap_samples

  # determine plot size and margins
  maxncharSamples <- max(sapply(bootstrap_samples,nchar))

  if(is.null(pwidth)) pwidth <- max(5,2+0.60*length(bootstrap_samples))
  if(is.null(pheight)) pheight <- 4+0.115*maxncharSamples
  if(is.null(mar)) {
    mar1 <- 0.65*maxncharSamples+1.2
    mar2 <- 4
    mar3 <- 3
    mar4 <- 2
    mar=c(mar1,mar2,mar3,mar4)
  }

  o <- order(hrdetect_res$hrdetect_output[,"Probability"],decreasing = TRUE)
  # jpeg(filename = paste0(outdir,"/HRDetect_bootstrap.jpg"),width = pwidth,height = pheight,res = pres)
  cairo_pdf(filename = paste0(outdir,"/HRDetect_bootstrap.pdf"),width = pwidth,height = pheight,pointsize = pointsize)
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
    points(myjitter, thisvalues[,o[i]], pch=20, col="#0067A51A")
    rect(xleft = i+0.15,ybottom = hrdetect_res$q_5_50_95[o[i],"5%"],xright = i+0.22,ytop = hrdetect_res$q_5_50_95[o[i],"95%"],col = "grey")
  }
  points(hrdetect_res$hrdetect_output[,"Probability"][o],pch=23,col="black",bg="white")

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
#' @param file_name name of the output file (pdf)
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

  maxncharSamples <- max(sapply(row.names(hrdetect_output),nchar,USE.NAMES = F))

  cairo_pdf(filename = file_name,
       width = 5+0.60*nrow(hrdetect_output),
       height = 6+0.115*maxncharSamples)
  par(mfrow=c(1,1),oma=c(0,0,1,0))
  mat <- matrix(c(1,2),ncol = 1)
  layout(mat, c(1), c(0.8,0.8+0.04*maxncharSamples))

  par(mar=c(1,7,3,15))
  bx <- barplot(height = hrdetect_output$Probability[reorder],
          ylim = c(0,1),ylab = "BRCA1/BRCA2\n def. score",
          names.arg = "",border = NA,
          cex.axis = 1.2,
          cex.lab=1.2,
          cex.names = 1.5,las=2)
  xpos <- as.vector(bx)
  if(length(xpos)>1){
    legendpos <- xpos[length(xpos)]*2 - xpos[length(xpos)-1]
    gap <- xpos[length(xpos)] - xpos[length(xpos)-1]
  }else{
    legendpos <- xpos[length(xpos)]*2
    gap <- xpos[length(xpos)]/2
  }

  legend_colors <- c('F38400', 'A1CAF1','875692','F3C300',   'BE0032', '8DB600')
  legend_colors <- paste0("#",legend_colors)
  legend(x=legendpos,y=1,legend = c("classification threshold"),lty = 2,col="red",cex = 1,bty = "n",xpd = TRUE)
  legend(x=legendpos+0.24,y=0.9,legend = c("deletion with MH",
                                 "Substitution Sig. 3",
                                 "Rearrangement Sig. 3",
                                 "Rearrangement Sig. 5",
                                 "HRD-LOH score",
                                 "Substitution Sig. 8"),fill =legend_colors,border = NA,
         cex = 1,bty = "n",xpd = TRUE)
  title(main=paste0(par_title), cex.main=1.3)
  # abline(a=0.7,b=0,col="red",lty = 2)
  lines(x=c(xpos[1]-gap/2,xpos[length(xpos)]+gap/2),y=c(0.7,0.7),col="red",lty = 2)
  # text(x = 0.75*nrow(hrdetect_output),y = 0.77,labels = "classification threshold",col = "red",cex = 1.5)
  matrix_to_plot <- t(hrdetect_output[reorder,col_needed])

  mar1 <- 0.65*maxncharSamples+1.2
  mar2 <- 7
  mar3 <- 1
  mar4 <- 15
  mar=c(mar1,mar2,mar3,mar4)
  # par(mar=c(10,7,1,15))
  par(mar=mar)
  barplot(height = matrix_to_plot,
          beside = TRUE,
          #names.arg = row.names(hrdetect_res),
          names.arg = rownames(hrdetect_output)[reorder],border = NA,
          las = 2,
          # legend.text = c("deletion with MH",
          #                 "Substitution Sig. 3",
          #                 "Rearrangement Sig. 3",
          #                 "Rearrangement Sig. 5",
          #                 "HRD-LOH score",
          #                 "Substitution Sig. 8"),
          # args.legend = list(x ='bottom', bty='n', inset=c(0,-1.0),horiz=TRUE,cex=0.9),
          # ylim = c(min(matrix_to_plot)*1.1,max(matrix_to_plot)*1.1),
          ylim = c(min(-1,min(matrix_to_plot)),max(1,max(matrix_to_plot))),
          # col=c("blue","red","black","green","orange","yellow"),
          col=legend_colors,
          ylab = "BRCA1/BRCA2\n def. contribution",
          cex.axis = 1.2,
          cex.lab = 1.2,
          cex.names = 1)

  dev.off()
}

