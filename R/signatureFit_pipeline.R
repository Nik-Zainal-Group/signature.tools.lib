# Andrea Degasperi ad923@cam.ac.uk, Serena Nik-Zainal group, University of Cambridge, UK, 2022

#' Signature fit pipeline
#'
#' This function is the main interface for computing signature fit using the signature.tools.lib R package.
#'
#' The pipeline will produce some feedback in the form or info, warning, and error messages.
#' Please check the output to see whether everything worked as planned.
#'
#' @param catalogues catalogues matrix, samples as columns, channels as rows. The mutation type of the catalogue will be inferred automatically by checking the rownames.
#' @param genome.v either "hg38" (will load BSgenome.Hsapiens.UCSC.hg38), "hg19" (will load BSgenome.Hsapiens.1000genomes.hs37d5), mm10 (will load BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10) or canFam3 (will load BSgenome.Cfamiliaris.UCSC.canFam3::BSgenome.Cfamiliaris.UCSC.canFam3)
#' @param organ If signatures is not specified, then use this paramenter to provide an organ name to automatically select appropriate signatures. Organ names and signature selection depends on the signature_version provided. When using RefSigv1 or RefSigv2 as signature_version organ-specific signatures will be used. Use one of the following organs: "Biliary", "Bladder", "Bone_SoftTissue", "Breast", "Cervix" (v1 only), "CNS", "Colorectal", "Esophagus", "Head_neck", "Kidney", "Liver", "Lung", "Lymphoid", "NET" (v2 only), "Oral_Oropharyngeal" (v2 only), "Ovary", "Pancreas", "Prostate", "Skin", "Stomach", "Uterus". If COSMICv2 or COSMICv3.2 are used, signatures are selected if the were found in the given organ/dataset. The mutation type is automatically inferred from the catalogue.
#' @param SNV_vcf_files list of file names corresponding to SNV VCF files to be used to construct 96-channel substitution catalogues. This should be a named vector, where the names indicate the sample name.
#' @param SNV_tab_files list of file names corresponding to SNV TAB files to be used to construct 96-channel substitution catalogues. This should be a named vector, where the names indicate the sample name. The files should contain a header in the first line with the following columns: chr, position, REF, ALT.
#' @param DNV_vcf_files list of file names corresponding to SNV/DNV VCF files to be used to construct 96-channel substitution catalogues. Adjacent SNVs will be combined into DNVs. This should be a named vector, where the names indicate the sample name.
#' @param DNV_tab_files list of file names corresponding to SNV/DNV TAB files to be used to construct 96-channel substitution catalogues. Adjacent SNVs will be combined into DNVs. This should be a named vector, where the names indicate the sample name. The files should contain a header in the first line with the following columns: chr, position, REF, ALT.
#' @param SV_bedpe_files list of file names corresponding to SV (Rearrangements) BEDPE files to be used to construct 32-channel rearrangement catalogues. This should be a named vector, where the names indicate the sample name. The files should contain a rearrangement for each row (two breakpoint positions should be on one row as determined by a pair of mates of paired-end sequencing) and should already be filtered according to the user preference, as all rearrangements in the file will be used and no filter will be applied. The files should contain a header in the first line with the following columns: "chrom1", "start1", "end1", "chrom2", "start2", "end2" and "sample" (sample name). In addition, either two columns indicating the strands of the mates, "strand1" (+ or -) and "strand2" (+ or -), or one column indicating the structural variant class, "svclass": translocation, inversion, deletion, tandem-duplication. The column "svclass" should correspond to (Sanger BRASS convention): inversion (strands +/- or -/+ and mates on the same chromosome), deletion (strands +/+ and mates on the same chromosome), tandem-duplication (strands -/- and mates on the same chromosome), translocation (mates are on different chromosomes).
#' @param signatures signatures should be a matrix or dataframe, signatures as columns, channels as rows. The mutation type of the signatures will be inferred automatically by checking the rownames. Use this parameter only if you want to use your own signatures. Leave NULL if you want to use the signatures provided by the package, for example by specifying a specific organ or signature_version.
#' @param rare_signatures used only when fit_method=FitMS, and the signature parameter is also given. The parameter rare_signatures should be a matrix or dataframe, signatures as columns, channels as rows. The mutation type of the signatures will be inferred automatically by checking the rownames.
#' @param signature_version either "COSMICv2", "COSMICv3.2", "RefSigv1" or "RefSigv2". If not specified, "RefSigv2 will be used. The mutation type is automatically inferred from the catalogue.
#' @param signature_names if no signatures have been provided using the signatures and rare_signatures parameters, and if no organ is specified, then signature_names can be used to specify a list of signature names, which should match the corresponding mutation type (inferred automatically) and reference signatures requested using the signature_version parameter.
#' @param fit_method either Fit or FitMS. Notice that automatic selection of signatures in FitMS is currrently available only for SNV mutations or catalogues, signature_version=RefSigv2 and specifying an organ. Alternatively, FitMS can be used by specifying both signatures (which will be considered common signatures) and rare_signatures parameters.
#' @param optimisation_method KLD or NNLS
#' @param useBootstrap set to TRUE to use bootstrap
#' @param nboot number of bootstraps to use, more bootstraps more accurate results
#' @param exposureFilterType use either fixedThreshold or giniScaledThreshold. When using fixedThreshold, exposures will be removed based on a fixed percentage with respect to the total number of mutations (threshold_percent will be used). When using giniScaledThreshold each signature will used a different threshold calculated as (1-Gini(signature))*giniThresholdScaling
#' @param threshold_percent threshold in percentage of total mutations in a sample, only exposures larger than threshold are considered
#' @param giniThresholdScaling scaling factor for the threshold type giniScaledThreshold, which is based on the Gini score of a signature
#' @param multiStepMode use one of the following: "constrainedFit", "partialNMF", "errorReduction", or "cossimIncrease".
#' @param residualNegativeProp maximum proportion of mutations (w.r.t. total mutations in a sample) that can be in the negative part of a residual when using the constrained least squares fit
#' when using multiStepMode=constrainedFit
#' @param minResidualMutations minimum number of mutations in a residual when using constrainedFit or partialNMF. Deactivated by default.
#' @param minCosSimRareSig minimum cosine similarity between a residual and a rare signature for considering the rare signature as a candidate for a sample when using constrainedFit or partialNMF
#' @param minErrorReductionPerc minimum percentage of error reduction for a signature to be considered as candidate when using the errorReduction method. The error is computed as mean absolute deviation
#' @param minCosSimIncrease minimum cosine similarity increase for a signature to be considered as candidate when using the cossimIncrease method
#' @param threshold_p.value p-value to determine whether an exposure is above the threshold_percent. In other words, this is the empirical probability that the exposure is lower than the threshold
#' @param rareSignatureTier either T1 or T2. For each organ we provide two lists of rare signatures that can be used. Tier 1 (T1) are rare signatures
#' that were observed in the requested organ. The problem with T1 is that it may be that a signature is not observed simply because there were not enough samples for a certain organ in the particular
#' dataset that was used to extract the signatures. So in general we advise to use Tier 2 (T2) signatures, which extend the rare signature to a wider number of rare signatures.
#' @param maxRareSigsPerSample masimum number of rare signatures that should be serched in each sample. In most situations, leaving this at 1 should be enough.
#' @param nparallel to use parallel specify >1
#' @param noFit if TRUE, terminate the pipeline early without running signature Fit. This is useful if one only wants to generate catalogues from mutation lists.
#' @param randomSeed set an integer random seed
#' @param verbose use FALSE to suppress messages
#' @return returns the fit object with activities/exposures of the signatures in the given sample and other information
#' @keywords mutational signatures fit
#' @export
#' @examples
#' res <- signatureFit_pipeline(catalogues,"Breast")
#' plotFitResults(res$fitResults,"results/")
signatureFit_pipeline <- function(catalogues=NULL,
                                  genome.v="hg19",
                                  organ=NULL,
                                  SNV_vcf_files=NULL,
                                  SNV_tab_files=NULL,
                                  DNV_vcf_files=NULL,
                                  DNV_tab_files=NULL,
                                  SV_bedpe_files=NULL,
                                  signatures=NULL,
                                  rare_signatures=NULL,
                                  signature_version="RefSigv2",
                                  signature_names=NULL,
                                  fit_method="FitMS", # either FitMS or Fit
                                  optimisation_method = "KLD",
                                  useBootstrap = FALSE,
                                  nboot = 200,
                                  exposureFilterType = "fixedThreshold", # or "giniScaledThreshold"
                                  threshold_percent = 5,
                                  giniThresholdScaling = 10,
                                  multiStepMode = "errorReduction", # or "partialNMF", or "errorReduction", or "cossimIncrease"
                                  threshold_p.value = 0.05,
                                  rareSignatureTier = "T2",  #either T1 for observed in organ or T2 for extended
                                  residualNegativeProp = 0.003,
                                  minResidualMutations = NULL,
                                  minCosSimRareSig = 0.8,
                                  minErrorReductionPerc = 15,
                                  minCosSimIncrease = 0.02,
                                  maxRareSigsPerSample = 1,
                                  noFit = FALSE,
                                  nparallel = 1,
                                  randomSeed = NULL,
                                  verbose = FALSE){

  message("[info signatureFit_pipeline] signatureFit pipeline starting!")

  #if multiple parallel cores are used, set it here
  doParallel::registerDoParallel(nparallel)

  if(!is.null(randomSeed)){
    set.seed(randomSeed)
  }

  # let's see if a catalogue was provided and whether we can infer the mutation type.
  mtype_catalogues <- NULL
  sampleNames <- NULL
  if(!is.null(catalogues)){
    mtype_catalogues <- getTypeOfMutationsFromChannels(catalogues)
    sampleNames <- colnames(catalogues)
  }

  # now check whether we need to build some catalogues
  specified_files <- ! c(is.null(SNV_vcf_files),
                         is.null(SNV_tab_files),
                        is.null(DNV_vcf_files),
                        is.null(DNV_tab_files),
                        is.null(SV_bedpe_files))
  names(specified_files) <- c("SNV_vcf_files",
                              "SNV_tab_files",
                              "DNV_vcf_files",
                              "DNV_tab_files",
                              "SV_bedpe_files")
  # In case the function that builds catalogues returns annotated mutations
  annotated_mutations <- NULL
  mtype_mutations <- NULL
  catalogues_mutations <- NULL

  if(sum(specified_files)>1){
    # too many mutation file types passed
    message("[error signatureFit_pipeline] too many mutation file types specified: ",paste(names(specified_files)[specified_files],collapse = ","),". Please specify only one.")
    return(NULL)
  }else if(sum(specified_files)==0){
    # If no mutation files were used, then just make sure the catalgues var is not NULL
    if(is.null(catalogues)){
      message("[error signatureFit_pipeline] no mutation files nor catalogues file specified. Nothing to run the analysis on.")
      return(NULL)
    }
  }else{
    # now here we have exaclty one mutation file to use, so build the catalogue accordingly
    if(!is.null(SNV_vcf_files)){
      message("[info signatureFit_pipeline] reading ",length(SNV_vcf_files)," SNV vcf mutation files and building catalogues.")

      cat_list <- foreach::foreach(i=1:length(SNV_vcf_files)) %dopar% {
        sample <- names(SNV_vcf_files)[i]
        res <- vcfToSNVcatalogue(SNV_vcf_files[sample],genome.v = genome.v)
        colnames(res$catalogue) <- sample
        res$muts <- cbind(data.frame(sampleName=rep(sample,nrow(res$muts)),stringsAsFactors = F),res$muts)
        res
      }
      catalogues_mutations <- data.frame(row.names = rownames(cat_list[[1]]$catalogue),stringsAsFactors = F)
      for(i in 1:length(cat_list)){
        catalogues_mutations <- cbind(catalogues_mutations,cat_list[[i]]$catalogue)
        annotated_mutations <- rbind(annotated_mutations,cat_list[[i]]$muts)
      }

    }else if(!is.null(SNV_tab_files)){
      message("[info signatureFit_pipeline] reading ",length(SNV_tab_files)," SNV tab mutation files and building catalogues.")

      cat_list <- foreach::foreach(i=1:length(SNV_tab_files)) %dopar% {
        sample <- names(SNV_tab_files)[i]
        subs <- read.table(file = SNV_tab_files[sample],
                           sep = "\t",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
        res <- tabToSNVcatalogue(subs,genome.v = genome.v)
        colnames(res$catalogue) <- sample
        res$muts <- cbind(data.frame(sampleName=rep(sample,nrow(res$muts)),stringsAsFactors = F),res$muts)
        res
      }
      catalogues_mutations <- data.frame(row.names = rownames(cat_list[[1]]$catalogue),stringsAsFactors = F)
      for(i in 1:length(cat_list)){
        catalogues_mutations <- cbind(catalogues_mutations,cat_list[[i]]$catalogue)
        annotated_mutations <- rbind(annotated_mutations,cat_list[[i]]$muts)
      }

    }else if(!is.null(DNV_vcf_files)){
      message("[info signatureFit_pipeline] reading ",length(DNV_vcf_files)," DNV vcf mutation files and building catalogues.")

      cat_list <- foreach::foreach(i=1:length(DNV_vcf_files)) %dopar% {
        sample <- names(DNV_vcf_files)[i]
        res <- vcfToDNVcatalogue(DNV_vcf_files[sample],genome.v = genome.v)
        colnames(res$DNV_catalogue) <- sample
        res$DNV_catalogue <- convertToAlexandrovChannels(res$DNV_catalogue)
        if(!is.null(res$DNV_table[["snv"]])) res$DNV_table[["snv"]] <- cbind(data.frame(sampleName=rep(sample,nrow(res$DNV_table[["snv"]])),stringsAsFactors = F),res$DNV_table[["snv"]])
        if(!is.null(res$DNV_table[["dnv"]])) res$DNV_table[["dnv"]] <- cbind(data.frame(sampleName=rep(sample,nrow(res$DNV_table[["dnv"]])),stringsAsFactors = F),res$DNV_table[["dnv"]])
        res
      }
      catalogues_mutations <- data.frame(row.names = rownames(cat_list[[1]]$DNV_catalogue),stringsAsFactors = F)
      annotated_mutations <- list()
      for(i in 1:length(cat_list)){
        catalogues_mutations <- cbind(catalogues_mutations,cat_list[[i]]$DNV_catalogue)
        if(!is.null(cat_list[[i]]$DNV_table[["snv"]])) annotated_mutations[["snv"]] <- rbind(annotated_mutations[["snv"]],cat_list[[i]]$DNV_table[["snv"]])
        if(!is.null(cat_list[[i]]$DNV_table[["dnv"]])) annotated_mutations[["dnv"]] <- rbind(annotated_mutations[["dnv"]],cat_list[[i]]$DNV_table[["dnv"]])
      }

    }else if(!is.null(DNV_tab_files)){
      message("[info signatureFit_pipeline] reading ",length(DNV_tab_files)," DNV tab mutation files and building catalogues.")

      cat_list <- foreach::foreach(i=1:length(DNV_tab_files)) %dopar% {
        sample <- names(DNV_tab_files)[i]
        muts <- read.table(file = DNV_tab_files[sample],
                           sep = "\t",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
        res <- tabToDNVcatalogue(muts)
        colnames(res$DNV_catalogue) <- sample
        res$DNV_catalogue <- convertToAlexandrovChannels(res$DNV_catalogue)
        if(!is.null(res$DNV_table[["snv"]])) res$DNV_table[["snv"]] <- cbind(data.frame(sampleName=rep(sample,nrow(res$DNV_table[["snv"]])),stringsAsFactors = F),res$DNV_table[["snv"]])
        if(!is.null(res$DNV_table[["dnv"]])) res$DNV_table[["dnv"]] <- cbind(data.frame(sampleName=rep(sample,nrow(res$DNV_table[["dnv"]])),stringsAsFactors = F),res$DNV_table[["dnv"]])
        res
      }
      catalogues_mutations <- data.frame(row.names = rownames(cat_list[[1]]$DNV_catalogue),stringsAsFactors = F)
      annotated_mutations <- list()
      for(i in 1:length(cat_list)){
        catalogues_mutations <- cbind(catalogues_mutations,cat_list[[i]]$DNV_catalogue)
        if(!is.null(cat_list[[i]]$DNV_table[["snv"]])) annotated_mutations[["snv"]] <- rbind(annotated_mutations[["snv"]],cat_list[[i]]$DNV_table[["snv"]])
        if(!is.null(cat_list[[i]]$DNV_table[["dnv"]])) annotated_mutations[["dnv"]] <- rbind(annotated_mutations[["dnv"]],cat_list[[i]]$DNV_table[["dnv"]])
      }

    }else if(!is.null(SV_bedpe_files)){
      message("[info signatureFit_pipeline] reading ",length(SV_bedpe_files)," SV bedpe mutation files and building catalogues.")

      cat_list <- foreach::foreach(i=1:length(SV_bedpe_files)) %dopar% {
        sample <- names(SV_bedpe_files)[i]
        sv_bedpe <- read.table(SV_bedpe_files[sample],sep = "\t",header = TRUE,
                               stringsAsFactors = FALSE,check.names = FALSE,comment.char = "")
        reslist <- bedpeToRearrCatalogue(sv_bedpe)
        # res <- reslist$rearr_catalogue
        # check that only one catalogue is generated. If not, take the one with more mutatations and raise a warning
        resncol <- ncol(reslist$rearr_catalogue)
        if(resncol>1){
          rescolsum <- apply(reslist$rearr_catalogue,2,sum)
          sampletouse <- which.max(rescolsum)[1]
          message(paste0("[warning signatureFit_pipeline] BEDPE file for sample ",sample," contained ",resncol," sample names: ",paste(colnames(res),collapse = ", "),". This could be due to germline rearrangements that should be removed. Using only the sample with the largest number of rearrangements (",colnames(res)[sampletouse],"). Please double check and rerun if necessary with only one sample for each BEDPE file."))
          reslist$rearr_catalogue <- reslist$rearr_catalogue[,sampletouse,drop=F]
          reslist$annotated_bedpe <- reslist$annotated_bedpe[reslist$annotated_bedpe$sample==colnames(reslist$rearr_catalogue)[sampletouse],,drop=F]
        }
        colnames(reslist$rearr_catalogue) <- sample
        reslist
      }
      catalogues_mutations <- data.frame(row.names = rownames(cat_list[[1]]$rearr_catalogue),stringsAsFactors = F)
      bedpecolumns <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2" , "sample","svclass","id", "is.clustered", "length")
      for(i in 1:length(cat_list)){
        catalogues_mutations <- cbind(catalogues_mutations,cat_list[[i]]$rearr_catalogue)
        annotated_mutations <- rbind(annotated_mutations,cat_list[[i]]$annotated_bedpe[,bedpecolumns,drop=F])
      }
    }

    # now we have built a catalogue table from the mutations
    mtype_mutations <- getTypeOfMutationsFromChannels(catalogues_mutations)
    if(is.null(catalogues)){
      catalogues <- catalogues_mutations
      mtype_catalogues <- mtype_mutations
      sampleNames <- colnames(catalogues)
    }else{
      # we already have a catalogue
      # let's check it is the same mutation type as the catalogue provided
      if(mtype_catalogues==mtype_mutations){
        # check for overlapping sample names and give precedence to the catalogue provided
        conflict_samples <- intersect(colnames(catalogues),colnames(catalogues_mutations))
        if(length(conflict_samples)>0){
          message("[warning signatureFit_pipeline] sample names conflict. The same sample names have been found in the catalogue provided and in the list of mutation files provided. Catalogues already provided will be used and mutations will be ignored for samples: ",paste(conflict_samples,collapse = ","),".")
          catalogues_mutations <- catalogues_mutations[,setdiff(colnames(catalogues_mutations),conflict_samples),drop=F]
        }
        # combine the catalogues:
        catalogues <- cbind(catalogues,catalogues_mutations)
        sampleNames <- colnames(catalogues)
      }else{
        message("[error signatureFit_pipeline] mutation files provided (",names(specified_files)[specified_files],
                ") has a mutation type (",mtype_mutations,") different from the mutation type of the catalogue (",mtype_catalogues,
                "). They need to be the same.")
        return(NULL)
      }
    }
  }

  # now we have the catalogues to fit and the annotated mutations ready
  # we should save these now in the return obj
  returnObj <- list()
  returnObj$catalogues <- catalogues
  returnObj$mtype_catalogues <- mtype_catalogues
  returnObj$annotated_mutations <- annotated_mutations

  # if no signature fit was requested then finish early and return catalogue and annotated mutations
  if(noFit){
    message("[info signatureFit_pipeline] Finishing pipeline early as no signature fit requested (noFit=TRUE). ",
            "Check the returned object for generated catalogues.")
    return(returnObj)
  }

  # find out what signatures you need, start from checking whether signatures were provided

  if(!is.null(signatures)){
    # check that the tumour type is the same as the catalogues'
    mtype_signatures <- getTypeOfMutationsFromChannels(signatures)
    if(mtype_catalogues!=mtype_signatures){
      message("[error signatureFit_pipeline] mutation files or catalogues provided have a mutation type (",
              mtype_catalogues,") different from the mutation type of the signatures file provided (",
              mtype_signatures,"). They need to be the same.")
      return(returnObj)
    }
  }

  if(!is.null(rare_signatures)){
    # check that the tumour type is the same as the catalogues'
    mtype_rare_signatures <- getTypeOfMutationsFromChannels(rare_signatures)
    if(mtype_catalogues!=mtype_rare_signatures){
      message("[error signatureFit_pipeline] mutation files or catalogues provided have a mutation type (",
              mtype_catalogues,") different from the mutation type of the rare signatures file provided (",
              mtype_rare_signatures,"). They need to be the same.")
      return(returnObj)
    }
  }

  # now, let's consider the following order/priorities:
  # 1. If the user provided signatures, then use those and ignore all other options
  #    If the user chose FitMS but did not provide rare signatures, then we should ask the user to provide the rare too
  # 2. If the user did not provide signatures, let's check if the user provided an organ.
  #    With an organ, we can select signatures automatically, depending on the signature_version parameter
  #    a) if signature_version is NULL, then assume we are using the latest organ-specific signatures
  #       (RefSigv2 if subs and dbs or RefSigv1 if rearr)
  #    b) if signature version is RefSigv1 or RefSigv2, then use organ specific signatures
  #    c) if signature version is COSMICv2 or COSMICv3.2, use a subset of signatures found in a given organ
  # 3. If no signatures are provided and no organ is specified, then use the COSMIC or RefSigs reference signatures.
  #    If a list of signatures is provided with signature_names, then use only the signatures in the list,
  #    while if signature_names is NULL then use all the signatures at once

  if(!is.null(signatures)){
    if(is.null(rare_signatures) & fit_method=="FitMS"){
      message("[error signatureFit_pipeline] user selected FitMS and provided a signatures file, which will be used for the common signatures. ",
              "Please also provide a file with the rare signatures using the parameter rare_signatures_file.")
      return(returnObj)
    }
    if(!is.null(organ)){
      message("[warning signatureFit_pipeline] organ parameter ",organ," will be ignored because the user has specified signature files.")
    }
    if(!is.null(signature_version)){
      message("[warning signatureFit_pipeline] signature_version parameter ",signature_version," will be ignored because the user has specified signature files.")
    }
  }else{
    # user has not provided signatures using file names, so let's see what was requested

    # let's check whether a specific signature_version was requested, and if not set to the latest
    if(is.null(signature_version)){
      message("[warning signatureFit_pipeline] signature_version parameter unspecified, using latest RefSigv2.")
      signature_version <- "RefSigv2"
    }

    # if FitMS was requested, but signature_version is not RefSigv2 or mtype_catalogues is not subs
    # or no organ is specified, then I can't use it, so revert to Fit
    if(fit_method=="FitMS"){
      if(signature_version!="RefSigv2" | mtype_catalogues!="subs" | is.null(organ)){
        message("[error signatureFit_pipeline] user selected FitMS and provided no signature files. ",
                "This means that common and rare signatures need to be automatically selected. ",
                "Automatic selection is currently available only for signature version RefSigv2, ",
                "substitutions mutation type or catalogues, and an organ should be specified. ",
                "Please correct your options or use fit_method=\"Fit\" instead.")
        return(returnObj)
      }
    }

    if(signature_version=="RefSigv2"){
      if(!is.null(organ)){
        # an organ was requested. If FitMS and subs were requested, then leave it to FitMS to get the appropriate signatures
        if(fit_method=="Fit"){
          if(mtype_catalogues %in% c("subs","DNV")){
            signatures <- getOrganSignatures(organ = organ,version = 2,typemut = mtype_catalogues)
          }else if(mtype_catalogues %in% c("rearr")){
            message("[warning signatureFit_pipeline] rearrangements RefSig mutational signatures only available in RefSigv1, ",
                    "switching to signature_version=RefSigv1")
            signatures <- getOrganSignatures(organ = organ,version = 1,typemut = mtype_catalogues)
          }else{
            message("[error signatureFit_pipeline] mutation type ",mtype_catalogues," not available for automatic selection of signatures. ",
                    "Please provide your own signatures using the signatures_file parameter, and also rare_signatures_file if using FitMS.")
            return(returnObj)
          }
          # need to check that the getOrganSignatures functions returned something
          if(ncol(signatures)==0){
            message("[error signatureFit_pipeline] RefSigv2 signatures not available for mutation type ",mtype_catalogues,
                    " and organ ",organ,".")
            return(returnObj)
          }
        }
      }else{
        # user did not provide an organ, use reference signatures and perhaps signature_names
        message("[info signatureFit_pipeline] no organ provided, so using the reference signatures appropriate for the signature version ",signature_version," and mutation type ",mtype_catalogues,".")
        if(mtype_catalogues == c("subs")){
          signatures <- referenceSignaturesSBSv2.03
        }else if(mtype_catalogues == c("DNV")){
          signatures <- referenceSignaturesDBSv1.01
        }else if(mtype_catalogues == c("rearr")){
          message("[warning signatureFit_pipeline] rearrangements RefSig mutational signatures only available in RefSigv1, ",
                  "switching to signature_version=RefSigv1")
          signatures <- RefSigv1_rearr
        }else{
          message("[error signatureFit_pipeline] mutation type ",mtype_catalogues," not available for automatic selection of signatures. ",
                  "Please provide your own signatures using the signatures_file parameter, and also rare_signatures_file if using FitMS.")
          return(returnObj)
        }
      }
    }else if(signature_version=="RefSigv1"){
      if(!is.null(organ)){
        # an organ was requested.
        if(fit_method=="Fit"){
          if(mtype_catalogues %in% c("subs","rearr")){
            signatures <- getOrganSignatures(organ = organ,version = 1,typemut = mtype_catalogues)
          }else if(mtype_catalogues %in% c("DNV")){
            message("[warning signatureFit_pipeline] DNV RefSig mutational signatures only available in RefSigv2, ",
                    "switching to signature_version=RefSigv2")
            signatures <- getOrganSignatures(organ = organ,version = 2,typemut = mtype_catalogues)
          }else{
            message("[error signatureFit_pipeline] mutation type ",mtype_catalogues," not available for automatic selection of signatures. ",
                    "Please provide your own signatures using the signatures_file parameter, and also rare_signatures_file if using FitMS.")
            return(returnObj)
          }
          # need to check that the getOrganSignatures functions returned something
          if(ncol(signatures)==0){
            message("[error signatureFit_pipeline] RefSigv1 signatures not available for mutation type ",mtype_catalogues,
                    " and organ ",organ,".")
            return(returnObj)
          }
        }
      }else{
        # user did not provide an organ, use reference signatures and perhaps signature_names
        message("[info signatureFit_pipeline] no organ provided, so using the reference signatures appropriate for the signature version ",signature_version," and mutation type ",mtype_catalogues,".")
        if(mtype_catalogues == c("subs")){
          signatures <- RefSigv1_subs
        }else if(mtype_catalogues == c("DNV")){
          message("[warning signatureFit_pipeline] DBS RefSig mutational signatures only available in RefSigv2, ",
                  "switching to signature_version=RefSigv2")
          signatures <- referenceSignaturesDBSv1.01
        }else if(mtype_catalogues == c("rearr")){
          signatures <- RefSigv1_rearr
        }else{
          message("[error signatureFit_pipeline] mutation type ",mtype_catalogues," not available for automatic selection of signatures. ",
                  "Please provide your own signatures using the signatures_file parameter, and also rare_signatures_file if using FitMS.")
          return(returnObj)
        }
      }
    }else if(signature_version=="COSMICv2"){
      if(!is.null(organ)){
        # TODO
        message("[error signatureFit_pipeline] organ signatures selection for COSMICv2 signatures not implemented yet. ",
                "Leave organ=NULL and select signatures manually with the signature_names parameter.")
        return(returnObj)

      }else{
        if(mtype_catalogues == c("subs")){
          signatures <- cosmic30
        }else{
          message("[error signatureFit_pipeline] using COSMICv2 signatures for ",mtype_catalogues," currenty not supported. ",
                  "Please provide your own signatures using the signatures_file parameter, and also rare_signatures_file if using FitMS.")
          return(returnObj)
        }

      }
    }else if(signature_version=="COSMICv3.2"){
      if(!is.null(organ)){
        # TODO
        message("[error signatureFit_pipeline] organ signatures selection for COSMICv3.2 signatures not implemented yet. ",
                "Leave organ=NULL and select signatures manually with the signature_names parameter.")
        return(returnObj)

      }else{
        if(mtype_catalogues == c("subs")){
          signatures <- COSMIC_v3.2_SBS_GRCh37
        }else{
          message("[error signatureFit_pipeline] using COSMICv2 signatures for ",mtype_catalogues," currenty not supported. ",
                  "Please provide your own signatures using the signatures_file parameter, and also rare_signatures_file if using FitMS.")
          return(returnObj)
        }
      }
    }else{
      message("[error signatureFit_pipeline] signature version ",signature_version," unrecognised.")
      return(returnObj)
    }
  }


  # signature selection using signature_names
  sigsToUseNames <- colnames(signatures)
  if(!is.null(signature_names)){
    if(fit_method=="Fit"){
      # check that the names make sense
      if(!all(signature_names %in% sigsToUseNames)){
        message("[error signatureFit_pipeline] using reference signatures for ",mtype_catalogues,", but some signature names provided ",
                "with the signature_names parameters do not seem to match the reference signatures. This is the list of ",
                "reference signatures available: ",paste(sigsToUseNames,collapse = ", "),".")
        return(returnObj)
      }
      # selected signatures only
      signatures <- signatures[,signature_names,drop=F]
    }else if(fit_method=="FitMS"){
      message("[warning signatureFit_pipeline] using FitMS, parameter signature_names will be ignored.")
    }

  }


  # Now, if fitMS was requested, and signature_type==NULL, then I think I can just call FitMS,
  # FitMS itself will take care of whether some of the parameters are invalid
  if(fit_method=="Fit"){
    message("[info signatureFit_pipeline] all set, running Fit.")
    fitRes <- Fit(catalogues = catalogues,
                  signatures = signatures,
                  exposureFilterType = exposureFilterType,
                  giniThresholdScaling = giniThresholdScaling,
                  threshold_percent = threshold_percent,
                  method = optimisation_method,
                  useBootstrap = useBootstrap,
                  nboot = nboot,
                  threshold_p.value = threshold_p.value,
                  nparallel = nparallel,
                  randomSeed = randomSeed,
                  verbose = verbose)
  }else if(fit_method=="FitMS"){
    message("[info signatureFit_pipeline] all set, running FitMS.")
    fitRes <- FitMS(catalogues = catalogues,
                    organ = organ,
                    rareSignatureTier = rareSignatureTier,
                    commonSignatures = signatures,
                    rareSignatures = rare_signatures,
                    method = optimisation_method,
                    exposureFilterType = exposureFilterType,
                    threshold_percent = threshold_percent,
                    giniThresholdScaling = giniThresholdScaling,
                    multiStepMode = multiStepMode,
                    residualNegativeProp = residualNegativeProp,
                    minResidualMutations = minResidualMutations,
                    minCosSimRareSig = minCosSimRareSig,
                    minErrorReductionPerc = minErrorReductionPerc,minCosSimIncrease = minCosSimIncrease,
                    useBootstrap = useBootstrap,
                    nboot = nboot,
                    threshold_p.value = threshold_p.value,
                    maxRareSigsPerSample = maxRareSigsPerSample,
                    nparallel = nparallel,
                    randomSeed = randomSeed,
                    verbose = verbose)
  }else{
    message("[error signatureFit_pipeline] unknown fit_method ",fit_method,".")
    return(returnObj)
  }

  message("[info signatureFit_pipeline] signatureFit pipeline completed!")

  # now just return the results
  returnObj$fitResults <- fitRes
  return(returnObj)

}

