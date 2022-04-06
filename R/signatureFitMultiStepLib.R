# Andrea Degasperi ad923@cam.ac.uk, Serena Nik-Zainal group, University of Cambridge, UK, 2021

#' Multi-Step Mutational Signatures Fit
#'
#' Given a set of mutational catalogues, this function will attempt fit mutational signature in a multi-step manner. In the first step, only the common signatures are fitted into the samples.
#' In the following steps, one or more rare signatures are fitted into the samples in addition to the common signatures. Common and rare signatures can be determined automatically by providing
#' the name of an organ, or can be supplied by the user.
#'
#' We provide four methods to identify the rare signatures in the samples:
#' "constrainedFit", "partialNMF", "errorReduction", or "cossimIncrease". The methods constrainedFit and partialNMF work in a similar way:
#' they identify a residual in each given sample, as the leftover mutations after fitting the common signatures. They will attempt to produce a mostly positive residual.
#' Each residual is then compared to each rare signature, and a rare signature is considered as a candidate rare signature for a sample if the cosine similarity between the residual
#' and the signature is at least minCosSimRareSig. One can also request that the residual is at least minResidualMutations. While constrainedFit will use a constrained least square fit where
#' the negative part of the residual is at most residualNegativeProp (a proportion of the number of mutations in the sample), partialNMF will instead use a few iterations of a KLD based NMF algorithm
#' where the matrix of the signatures contains the common signature and an additional signatures that needs to be estimated (NNLM package). The methods errorReduction and cossimIncrease work
#' in a similar way: they will fit the common signatures along with one additional rare signature, testing all rare signatures one at a time, and then determine difference in error
#' (or cosine similarity) between fitting the common signatures only and with the additional rare signatures. If the error reduction is at least minErrorReductionPerc (or the cosine similarity
#' increase is at least minCosSimIncrease then the rare signature will be considered as a candidate.
#'
#' After any of ghe procedures above, each sample may have multiple candidate rare signatures, so one is chosen according to the highest associated cosine similarity either of
#' the residual to the candidate rare signature (constrainedFit and partialNMF methods), or of the catalogue and the reconstructed sample (errorReduction and cossimIncrease methods).
#' It is then possible to plot all the fits with plotFitMS and even change the choise of the candidate rare signature using the function fitMerge.
#'
#' A post fit exposure filter will reduce the false positive singature assignments by setting to zero exposure values that
#' are below a certain threshold. We provide two exposureFilterType methods: fixedThreshold and giniScaledThreshold. The
#' fixedThreshold method will set to zero exposures that are below a fixed threshold given as a percentage of the mutations
#' in a sample (parameter threshold_percent), while the method giniScaledThreshold will use a different threshold for each
#' signature, computed as (1-Gini(signature))*giniThresholdScaling, which will also be a percentage of the mutations in a sample.
#'
#' @param catalogues catalogues matrix, samples as columns, channels as rows
#' @param organ #automatically sets the commonSignatures and rareSignatures parameters, which can be left as NULL. The following organs are available:
#' "Biliary", "Bladder", "Bone_SoftTissue", "Breast", "CNS", "Colorectal", "Esophagus", "Head_neck", "Kidney", "Liver", "Lung", "Lymphoid", "Myeloid",
#' "NET", "Oral_Oropharyngeal", "Ovary", "Pancreas", "Prostate", "Skin", "Stomach", "Uterus"
#' @param rareSignatureTier either T1 or T2. For each organ we provide two lists of rare signatures that can be used. Tier 1 (T1) are rare signatures
#' that were observed in the requested organ. The problem with T1 is that it may be that a signature is not observed simply because there were not enough samples for a certain organ in the particular
#' dataset that was used to extract the signatures. So in general we advise to use Tier 2 (T2) signatures, which extend the rare signature to a wider number of rare signatures.
#' @param commonSignatures signatures, signatures as columns, channels as rows. These are the signatures that are assumed to be present in most samples and will be used in the first step.
#' Can be set automatically by specifying the organ parameter
#' @param rareSignatures signatures, signatures as columns, channels as rows. These are the signatures that are assumed to be rarely present in a sample, at most maxRareSigsPerSample rare signatures in each sample.
#' Can be set automatically by specifying the organ parameter and the rareSignatureTier parameter
#' @param method KLD or NNLS
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
#' @param useBootstrap set to TRUE to use bootstrap
#' @param nboot number of bootstraps to use, more bootstraps more accurate results
#' @param threshold_p.value p-value to determine whether an exposure is above the threshold_percent. In other words, this is the empirical probability that the exposure is lower than the threshold
#' @param maxRareSigsPerSample masimum number of rare signatures that should be serched in each sample. In most situations, leaving this at 1 should be enough.
#' @param nparallel to use parallel specify >1
#' @param randomSeed set an integer random seed
#' @param verbose use FALSE to suppress messages
#' @return returns the activities/exposures of the signatures in the given sample and other information
#' @keywords mutational signatures fit
#' @export
#' @references A. Degasperi, X. Zou, T. D. Amarante, ..., H. Davies, Genomics England Research Consortium, S. Nik-Zainal. Substitution mutational signatures in whole-genome-sequenced cancers in the UK population. Science, 2022.
#' @examples
#' res <- FitMS(catalogues,"Breast")
#' plotFitMS(res,"results/")
FitMS <- function(catalogues,
                  organ = NULL, #automatically sets the common and rare signatures
                  rareSignatureTier = "T2",  #either T1 for observed in organ or T2 for extended
                  commonSignatures = NULL,
                  rareSignatures = NULL,
                  method = "KLD",
                  exposureFilterType = "fixedThreshold", # or "giniScaledThreshold"
                  threshold_percent = 5,
                  giniThresholdScaling = 10,
                  multiStepMode = "errorReduction", # or "partialNMF", or "errorReduction", or "cossimIncrease"
                  residualNegativeProp = 0.003,
                  minResidualMutations = NULL,
                  minCosSimRareSig = 0.8,
                  minErrorReductionPerc = 15,
                  minCosSimIncrease = 0.02,
                  useBootstrap = FALSE,
                  nboot = 200,
                  threshold_p.value = 0.05,
                  maxRareSigsPerSample = 1,
                  nparallel = 1,
                  randomSeed = NULL,
                  verbose = FALSE){

  # check type of mutations
  typeofmuts <- getTypeOfMutationsFromChannels(catalogues)

  # check common signatures
  if(!is.null(commonSignatures)){
    if(verbose) message("[info FitMS] Using user provided commonSignatures.")
  }else{
    if(is.null(organ)){
      message("[error FitMS] No organ was specified and no commonSignatures were given. Nothing to do.")
      return(NULL)
    }else{
      if(typeofmuts=="subs"){
        commonSignatures <- organSignaturesSBSv2.03[,strsplit(sigsForFittingSBSv2.03[organ,"common"],split = ",")[[1]],drop=F]
        if(ncol(commonSignatures)==0) {
          message("[error FitMS] No common signatures is associated with organ ",organ,". Nothing to do.")
          return(NULL)
        }
      }else{
        message("[error FitMS] Type of mutations ",typeofmuts," not available, please provide your own commonSignatures. Nothing to do.")
        return(NULL)
      }
    }
  }

  # check rare signatures
  if(!is.null(rareSignatures)){
    if(verbose) message("[info FitMS] Using user provided rareSignatures.")
  }else{
    if(rareSignatureTier %in% c("T1","T2")){
      if(is.null(organ)){
        message("[error FitMS] No organ was specified and no rareSignatures were given. Nothing to do.")
        return(NULL)
      }else{
        if(typeofmuts=="subs"){
          rareSignatures <- referenceSignaturesSBSv2.03[,strsplit(sigsForFittingSBSv2.03[organ,paste0("rare",rareSignatureTier)],split = ",")[[1]],drop=F]
          if(ncol(rareSignatures)==0) {
            message("[error FitMS] No rare signatures are associated with organ ",organ,". Nothing to do.")
            return(NULL)
          }
        }else{
          message("[error FitMS] Type of mutations ",typeofmuts," not available, please provide your own rareSignatures. Nothing to do.")
          return(NULL)
        }
      }
    }else{
      message("[error FitMS] invalid rareSignatureTier ",rareSignatureTier,". Please provide your own rareSignatures or choose between T1 and T2.")
      return(NULL)
    }
  }

  # ---- now determine which samples are likely to have rare signatures and how many
  # after this section, all we need are:
  # - whichSamplesMayHaveRareSigs
  # - candidateRareSigs
  # - candidateRareSigsCosSim

  commonSigsOnlyCosSim <- list()
  commonSigsOnlyError <- list()
  whichSamplesMayHaveRareSigs <- c()
  candidateRareSigs <- list()
  candidateRareSigsCosSim <- list()
  candidateRareSigsError <- list()

  # if you set maxRareSigsPerSample to 0, you should get the common sig fits only
  if(maxRareSigsPerSample > 0){

    for (depth in 1:maxRareSigsPerSample) {
      # depth <- 1

      for(i in 1:ncol(catalogues)){
        # i <- 1
        # if this is depth==1 then I am going to test all samples
        # if this is depth>1 then I am going to test only the samples that have depth-1 rare signatures
        # i.e. I am just going to search for one additional signature

        sampleName <- colnames(catalogues)[i]
        currentCatalogue <- catalogues[,sampleName,drop=F]
        ncandidatesAtPreviousDepth <- 0
        candidatesAtPreviousDepth <- NULL
        if(depth>1) {
          # check if something was found at previous depth for this sample
          candidatesAtPreviousDepth <- unlist(sapply(candidateRareSigs[[sampleName]],function(x){
            if(length(strsplit(x,split = ":")[[1]])==(depth-1)){
              return(x)
            }else{
              return(NULL)
            }
          },USE.NAMES = F))
          if(!is.null(candidatesAtPreviousDepth)) ncandidatesAtPreviousDepth <- length(candidatesAtPreviousDepth)
        }

        if((depth==1 | ncandidatesAtPreviousDepth>0) & (ncol(rareSignatures)-depth>=0)){
          # OK we are now GO for this sample to search for additional signatures
          # we are searching either only once if depth is 1, or as many times as necessary if depth is >1
          # also we know we have enough rare signatures to search up to this depth
          nloop <- 1
          if(depth>1) nloop <- ncandidatesAtPreviousDepth
          for(loopi in 1:nloop){
            # loopi <- 1
            if(depth==1){
              currentRareSigs <- NULL
              commonSigsToUse <- commonSignatures
              rareSigsToUse <- rareSignatures
            }else if(depth>1){
              currentRareSigs <- strsplit(candidatesAtPreviousDepth[loopi],split = ":")[[1]]
              # add the current rare signatures to the common set and remove them from the rare set
              commonSigsToUse <- cbind(commonSignatures,rareSignatures[,currentRareSigs,drop=F])
              rareSigsToUse <- rareSignatures[,setdiff(colnames(rareSignatures),currentRareSigs),drop=F]
            }

            # now we search with the various methods
            if(multiStepMode=="partialNMF" | multiStepMode=="constrainedFit"){
              #
              # compute residuals
              E <- flexconstr_sigfit_multipleSamples(as.matrix(commonSigsToUse),as.matrix(currentCatalogue),allmut_tolratio = residualNegativeProp)
              if(multiStepMode=="constrainedFit"){
                R <- currentCatalogue - as.matrix(commonSigsToUse) %*% as.matrix(E)
              }else if(multiStepMode=="partialNMF"){
                R <- partialNMFresidual(currentCatalogue,commonSigsToUse,E,max.iter = 50)
              }

              # check for min residual
              addSig <- TRUE
              if(!is.null(minResidualMutations)) addSig <- sum(R)>=minResidualMutations

              # now get which rare signatures are the most similar
              positiveR <- R
              positiveR[positiveR<0] <- 0
              simMatrix <- computeCorrelationOfTwoSetsOfSigs(positiveR,rareSigsToUse)

              # there may be more than one signature suitable
              simMatrix <- unlist(simMatrix)
              pos <- which(simMatrix>=minCosSimRareSig)
              if(length(pos)!=0 & addSig){
                if(depth==1){
                  candidateRareSigs[[colnames(catalogues)[i]]] <- colnames(rareSigsToUse)[pos]
                  candidateRareSigsCosSim[[colnames(catalogues)[i]]] <- simMatrix[pos]
                  names(candidateRareSigsCosSim[[colnames(catalogues)[i]]]) <- candidateRareSigs[[colnames(catalogues)[i]]]
                  whichSamplesMayHaveRareSigs <- c(whichSamplesMayHaveRareSigs,i)
                }else if(depth>1){
                  candidateRareSigs[[colnames(catalogues)[i]]] <- c(candidateRareSigs[[colnames(catalogues)[i]]],paste(currentRareSigs,colnames(rareSigsToUse)[pos],sep=":"))
                  candidateRareSigsCosSim[[colnames(catalogues)[i]]] <- c(candidateRareSigsCosSim[[colnames(catalogues)[i]]],simMatrix[pos])
                  names(candidateRareSigsCosSim[[colnames(catalogues)[i]]]) <- candidateRareSigs[[colnames(catalogues)[i]]]
                }
              }

            }else if(multiStepMode=="errorReduction" | multiStepMode=="cossimIncrease"){

              # get first fit with all
              quickFitCommon <- Fit(catalogues = currentCatalogue,
                                      signatures = commonSigsToUse,
                                      nboot = nboot,
                                      exposureFilterType = "fixedThreshold",
                                      giniThresholdScaling = giniThresholdScaling,
                                      threshold_percent = -1,
                                      threshold_p.value = threshold_p.value,
                                      method = method,
                                      useBootstrap = FALSE,
                                      nparallel = nparallel)

              quickFit <- list()
              for(j in 1:ncol(rareSigsToUse)){
                # j <- 1
                if(verbose) message("[info FitMS] multiStepMode ",multiStepMode,": depth ",depth,", sample ",i," of ",ncol(catalogues)," fitting rare signature ",j," of ",ncol(rareSigsToUse),"")
                quickFit[[j]] <- Fit(catalogues = currentCatalogue,
                                       signatures = cbind(commonSigsToUse,rareSigsToUse[,j,drop=F]),
                                       nboot = nboot,
                                       exposureFilterType = "fixedThreshold",
                                       giniThresholdScaling = giniThresholdScaling,
                                       threshold_percent = -1,
                                       threshold_p.value = threshold_p.value,
                                       method = method,
                                       useBootstrap = FALSE,
                                       nparallel = nparallel)
              }

              # commonError <- KLD(currentCatalogue,as.matrix(commonSigsToUse) %*% as.matrix(quickFitCommon$exposures))
              commonError <- MAD(currentCatalogue,as.matrix(commonSigsToUse) %*% as.matrix(t(quickFitCommon$exposures[,1:(ncol(quickFitCommon$exposures)-1),drop=F])))
              commonCosSim <- cos_sim(currentCatalogue,as.matrix(commonSigsToUse) %*% as.matrix(t(quickFitCommon$exposures[,1:(ncol(quickFitCommon$exposures)-1),drop=F])))
              if(depth==1){
                commonSigsOnlyCosSim[[colnames(catalogues)[i]]] <- commonCosSim
                commonSigsOnlyError[[colnames(catalogues)[i]]] <- commonError
              }

              allError <- sapply(1:length(quickFit), function(j) {
                MAD(currentCatalogue,as.matrix(cbind(commonSigsToUse,rareSigsToUse[,j,drop=F])) %*% as.matrix(t(quickFit[[j]]$exposures[,1:(ncol(quickFit[[j]]$exposures)-1),drop=F])))
              })
              allCosSim <- sapply(1:length(quickFit), function(j) {
                cos_sim(currentCatalogue,as.matrix(cbind(commonSigsToUse,rareSigsToUse[,j,drop=F])) %*% as.matrix(t(quickFit[[j]]$exposures[,1:(ncol(quickFit[[j]]$exposures)-1),drop=F])))
              })
              # error
              errorRed <- (commonError-allError)/commonError*100
              cossimIncr <- allCosSim-commonCosSim
              if(multiStepMode=="errorReduction"){
                pos <- which(errorRed>=minErrorReductionPerc)
              }else if(multiStepMode=="cossimIncrease"){
                pos <- which(cossimIncr>=minCosSimIncrease)
              }

              if(length(pos)!=0){
                if(depth==1){
                  candidateRareSigs[[colnames(catalogues)[i]]] <- colnames(rareSigsToUse)[pos]
                  candidateRareSigsCosSim[[colnames(catalogues)[i]]] <- allCosSim[pos]
                  candidateRareSigsError[[colnames(catalogues)[i]]] <- allError[pos]
                  names(candidateRareSigsCosSim[[colnames(catalogues)[i]]]) <- candidateRareSigs[[colnames(catalogues)[i]]]
                  names(candidateRareSigsError[[colnames(catalogues)[i]]]) <- candidateRareSigs[[colnames(catalogues)[i]]]
                  whichSamplesMayHaveRareSigs <- c(whichSamplesMayHaveRareSigs,i)
                }else if(depth>1){
                  candidateRareSigs[[colnames(catalogues)[i]]] <- c(candidateRareSigs[[colnames(catalogues)[i]]],paste(currentRareSigs,colnames(rareSigsToUse)[pos],sep=":"))
                  candidateRareSigsCosSim[[colnames(catalogues)[i]]] <- c(candidateRareSigsCosSim[[colnames(catalogues)[i]]],allCosSim[pos])
                  candidateRareSigsError[[colnames(catalogues)[i]]] <- c(candidateRareSigsError[[colnames(catalogues)[i]]],allError[pos])
                  names(candidateRareSigsCosSim[[colnames(catalogues)[i]]]) <- candidateRareSigs[[colnames(catalogues)[i]]]
                  names(candidateRareSigsError[[colnames(catalogues)[i]]]) <- candidateRareSigs[[colnames(catalogues)[i]]]
                }
              }

            }

            # --- done with the methods
          }
        }

      }

    }
  }

  # at this point it would be good to have a remove duplicates step
  allcomb <- function(x){
    tmpx <- strsplit(x,split = ":")[[1]]
    sapply(combinat::permn(tmpx),function(y) paste(y,collapse = ":"))
  }

  tmpcandidateRareSigs <- list()
  tmpcandidateRareSigsCosSim <- list()
  tmpcandidateRareSigsError <- list()
  for (s in names(candidateRareSigs)){
    # s <- names(candidateRareSigs)[1]
    for (si in 1:length(candidateRareSigs[[s]])){
      # si <- 1
      if(!any(allcomb(candidateRareSigs[[s]][si]) %in% tmpcandidateRareSigs[[s]])){
        tmpcandidateRareSigs[[s]] <- c(tmpcandidateRareSigs[[s]],candidateRareSigs[[s]][si])
        tmpcandidateRareSigsCosSim[[s]] <- c(tmpcandidateRareSigsCosSim[[s]],candidateRareSigsCosSim[[s]][si])
        tmpcandidateRareSigsError[[s]] <- c(tmpcandidateRareSigsError[[s]],candidateRareSigsError[[s]][si])
      }
    }
  }
  candidateRareSigs <- tmpcandidateRareSigs
  candidateRareSigsCosSim <- tmpcandidateRareSigsCosSim
  candidateRareSigsError <- tmpcandidateRareSigsError

  samples <- list()

  for(i in 1:ncol(catalogues)){
    # i <- 1
    sampleName <- colnames(catalogues)[i]
    currentCatalogue <- catalogues[,sampleName,drop=F]
    samples[[sampleName]] <- list()
    samples[[sampleName]]$catalogue <- currentCatalogue

    if(verbose) message("[info FitMS] fitting sample ",i," of ",ncol(catalogues),": ",sampleName)

    # begin by fitting
    samples[[sampleName]]$fitCommonOnly <- Fit(catalogues = currentCatalogue,
                                                 signatures = commonSignatures,
                                                 nboot = nboot,
                                                 exposureFilterType = exposureFilterType,
                                                 giniThresholdScaling = giniThresholdScaling,
                                                 threshold_percent = threshold_percent,
                                                 threshold_p.value = threshold_p.value,
                                                 method = method,
                                                 useBootstrap = useBootstrap,
                                                 nparallel = nparallel,
                                                 randomSeed = randomSeed)

    # now check if I should try to fit some rare signatures as well
    if(i %in% whichSamplesMayHaveRareSigs){
      samples[[sampleName]]$fitWithRare <- list()

      for (j in 1:length(candidateRareSigs[[sampleName]])){
        # j <- 1
        # consider the case in which >1 signatures are possible and separated by :
        currentRareSigName <- strsplit(candidateRareSigs[[sampleName]][j],split = ":")[[1]]
        currentRareSig <- rareSignatures[,currentRareSigName,drop=F]
        samples[[sampleName]]$fitWithRare[[paste(currentRareSigName,collapse = ":")]] <- Fit(catalogues = currentCatalogue,
                                                                                               signatures = cbind(commonSignatures,currentRareSig),
                                                                                               nboot = nboot,
                                                                                               exposureFilterType = exposureFilterType,
                                                                                               giniThresholdScaling = giniThresholdScaling,
                                                                                               threshold_percent = threshold_percent,
                                                                                               threshold_p.value = threshold_p.value,
                                                                                               method = method,
                                                                                               useBootstrap = useBootstrap,
                                                                                               nparallel = nparallel,
                                                                                               randomSeed = randomSeed)
      }
    }
  }

  # collect all info
  resObj <- list()
  # add data and options used
  resObj$fitAlgorithm <- "FitMS"
  resObj$catalogues <- catalogues
  resObj$organ <- organ
  resObj$commonSignatures <- commonSignatures
  resObj$rareSignatures <- rareSignatures
  resObj$method <- method
  resObj$exposureFilterType <- exposureFilterType
  resObj$multiStepMode <- multiStepMode
  resObj$giniThresholdScaling <- giniThresholdScaling
  resObj$minErrorReductionPerc <- minErrorReductionPerc
  resObj$minCosSimIncrease <-minCosSimIncrease
  resObj$threshold_percent <- threshold_percent
  resObj$residualNegativeProp <- residualNegativeProp
  resObj$minResidualMutations <- minResidualMutations
  resObj$minCosSimRareSig <- minCosSimRareSig
  resObj$useBootstrap <- useBootstrap
  resObj$nboot <- nboot
  resObj$threshold_p.value <- threshold_p.value
  resObj$rareSignatureTier <- rareSignatureTier
  resObj$maxRareSigsPerSample <- maxRareSigsPerSample
  # add all fit results and info
  resObj$whichSamplesMayHaveRareSigs <- colnames(catalogues)[whichSamplesMayHaveRareSigs]
  resObj$candidateRareSigs <- candidateRareSigs
  resObj$candidateRareSigsCosSim <- candidateRareSigsCosSim
  resObj$candidateRareSigsError <- candidateRareSigsError
  resObj$commonSigsOnlyError <- commonSigsOnlyError
  resObj$commonSigsOnlyCosSim <- commonSigsOnlyCosSim
  resObj$samples <- samples
  # compute fit merge
  resObj <- fitMerge(resObj)
  return(resObj)
}


#' Signature Fit
#'
#' This function provides basic signature fit functionalities.
#' Fit a given set of mutational signatures into mutational catalogues to estimate
#' the activty/exposure of each of the given signatures in the catalogues.
#'
#' This is a standard interface to signature fit functions with/without bootstrap. The object returned by this
#' function can be passed to the plotFit() function for automated plotting of the results.
#'
#' A post fit exposure filter will reduce the false positive singature assignments by setting to zero exposure values that
#' are below a certain threshold. We provide two exposureFilterType methods: fixedThreshold and giniScaledThreshold. The
#' fixedThreshold method will set to zero exposures that are below a fixed threshold given as a percentage of the mutations
#' in a sample (parameter threshold_percent), while the method giniScaledThreshold will use a different threshold for each
#' signature, computed as (1-Gini(signature))*giniThresholdScaling, which will also be a percentage of the mutations in a sample.
#'
#' @param catalogues catalogues matrix, samples as columns, channels as rows
#' @param signatures mutational signatures to bw fitted into the sample catalgues, signatures as columns and channels as rows
#' @param method KLD or NNLS
#' @param exposureFilterType use either fixedThreshold or giniScaledThreshold. When using fixedThreshold, exposures will be removed based on a fixed percentage with respect to the total number of mutations (threshold_percent will be used). When using giniScaledThreshold each signature will used a different threshold calculated as (1-Gini(signature))*giniThresholdScaling
#' @param threshold_percent threshold in percentage of total mutations in a sample, only exposures larger than threshold are considered
#' @param giniThresholdScaling scaling factor for the threshold type giniScaledThreshold, which is based on the Gini score of a signature
#' @param useBootstrap set to TRUE to use the signature fit with bootstrap method
#' @param nboot number of bootstraps to use, more bootstraps more accurate results (use only when useBootstrap=TRUE)
#' @param threshold_p.value p-value to determine whether an exposure is above the threshold_percent.
#' In other words, this is the empirical probability that the exposure is lower than the threshold (use only when useBootstrap=TRUE)
#' @param nparallel to use parallel specify >1
#' @param randomSeed set an integer random seed (use only when useBootstrap=TRUE)
#' @param verbose use FALSE to suppress messages
#' @return returns the activities/exposures of the signatures in the given sample and other information
#' @keywords mutational signatures fit
#' @export
#' @examples
#' res <- Fit(catalogues,getOrganSignatures("Breast"))
#' plotFit(res,"results/")
Fit <- function(catalogues,
                signatures,
                exposureFilterType = "fixedThreshold", # or "giniScaledThreshold"
                giniThresholdScaling = 10,
                threshold_percent = 5,
                method = "KLD",
                useBootstrap = FALSE,
                nboot = 200,
                threshold_p.value = 0.05,
                nparallel = 1,
                randomSeed = NULL,
                verbose = FALSE){
  # initialise return object
  fitRes <- list()
  fitRes$fitAlgorithm <- "Fit"

  if(useBootstrap){

    # compute the full results
    resFitBoot <- SignatureFit_withBootstrap(cat = catalogues,
                                             signature_data_matrix = signatures,
                                             nboot = nboot,
                                             exposureFilterType = exposureFilterType,
                                             threshold_percent = threshold_percent,
                                             giniThresholdScaling = giniThresholdScaling,
                                             threshold_p.value = threshold_p.value,
                                             method = method,
                                             verbose = verbose,
                                             doRound = FALSE,
                                             nparallel = nparallel,
                                             randomSeed = randomSeed,
                                             showDeprecated = F)

    fitRes$exposures <- resFitBoot$E_median_filtered
    fitRes$unassigned_muts <- apply(catalogues,2,sum) - apply(fitRes$exposures,2,sum)
    fitRes$unassigned_muts_perc <- fitRes$unassigned_muts/apply(catalogues,2,sum)*100
    fitRes$bootstrap_exposures <- resFitBoot$boot_list
    fitRes$bootstrap_exposures_samples <- resFitBoot$samples_list
    fitRes$bootstrap_exposures_pvalues <- resFitBoot$E_p.values
    fitRes$nboot <- nboot
    fitRes$threshold_p.value <- threshold_p.value

  }else{
    resFit <- SignatureFit(cat =catalogues,
                           signature_data_matrix = signatures,
                           method = method,
                           doRound = FALSE,
                           verbose = verbose,
                           showDeprecated = F)
    if(exposureFilterType=="giniScaledThreshold"){
      sigInvGini <- 1 - apply(signatures,2,giniCoeff)
      giniThresholdPerc <- giniThresholdScaling*sigInvGini
      # set to zero differently for each signature
      for(i in 1:length(giniThresholdPerc)) resFit[i,resFit[i,]/apply(catalogues,2,sum)*100<giniThresholdPerc[i]] <- 0
    }else if(exposureFilterType=="fixedThreshold"){
      resFit[resFit/matrix(apply(catalogues,2,sum),ncol = ncol(resFit),nrow = nrow(resFit),byrow = T)*100<threshold_percent] <- 0
    }

    fitRes$exposures <- resFit
    fitRes$unassigned_muts <- apply(catalogues,2,sum) - apply(fitRes$exposures,2,sum)
    fitRes$unassigned_muts_perc <- fitRes$unassigned_muts/apply(catalogues,2,sum)*100
    fitRes$bootstrap_exposures <- NA
    fitRes$bootstrap_exposures_samples <- NA
    fitRes$bootstrap_exposures_pvalues <- NA
    fitRes$nboot <- NA
    fitRes$threshold_p.value <- NA

  }

  fitRes$useBootstrap <- useBootstrap
  fitRes$exposureFilterType <- exposureFilterType
  fitRes$giniThresholdScaling <- ifelse(exposureFilterType=="giniScaledThreshold",giniThresholdScaling,NA)
  fitRes$threshold_percent <- ifelse(exposureFilterType=="fixedThreshold",threshold_percent,NA)
  fitRes$catalogues <- catalogues
  fitRes$signatures <- signatures
  fitRes$method <- method
  fitRes$exposures <- t(rbind(fitRes$exposures,unassigned=fitRes$unassigned_muts))

  return(fitRes)
}


#' Selecting and merging signature fit results from Multi-Step signature fit
#'
#' This function is used to select candidate rare signature solutions from a result object obtained
#' using the multi-step signature fit FitMS function.
#'
#' When running FitMS, some samples may have multiple candidate rare signatures or rare signature combinations that
#' fit the sample reasonably well. This function selects the best rare signature for each sample based on cosine similarity and returns
#' an updated object with a summary exposures composed by each selected fit solution for each sample.
#'
#' The choise of rare signature can be changed using the parameter forceRareSigChoice.
#'
#' @param resObj result object obtained from the FitMS function
#' @param forceRareSigChoice if NULL this function will select the rare signature candidate with the highest associated cosine similarity.
#' If no rare signature is found, then the solution with only the common signatures is selected.
#' To select specific candidates, specify them in the forceRareSigChoice list object, in the form forceRareSigChoice[["sample_name"]] <- "rareSigName".
#' To select the solution with only the common signatures for a sample use forceRareSigChoice[["sample_name"]] <- "common"
#' @return returns the updated resObj object with updated exposures and rareSigChoice objects.
#' If bootstrap was used, bootstraps of selected solutions can be found in the variable bootstrap_exposures_samples
#' @export
fitMerge <- function(resObj,forceRareSigChoice=NULL){
  # build exposure matrix with all signatures in the columns
  # will remove the signatures not present at the end
  rareSigChoice <- list()
  exposures_merge <- matrix(0,ncol = ncol(resObj$commonSignatures)+ncol(resObj$rareSignatures)+1,nrow = ncol(resObj$catalogues),
                            dimnames = list(colnames(resObj$catalogues),c(colnames(resObj$commonSignatures),colnames(resObj$rareSignatures),"unassigned")))
  bootstrap_exposures_samples <- list()
  for (i in 1:ncol(resObj$catalogues)){
    # i <- 1
    currentSample <- colnames(resObj$catalogues)[i]

    # check for attempts to force rare signature choice
    forceRareSig <- NULL
    forceCommon <- FALSE
    if(!is.null(forceRareSigChoice[[currentSample]])){
      if(forceRareSigChoice[[currentSample]]=="common"){
        forceCommon <- TRUE
      }else if(!is.null(resObj$samples[[currentSample]]$fitWithRare[[forceRareSigChoice[[currentSample]]]])){
        forceRareSig <- forceRareSigChoice[[currentSample]]
      }else{
        message("[warning fitMerge] Attempt to force rare sig ",forceRareSigChoice[[currentSample]]," for sample ",currentSample," failed because the rare signature is not in the fitWithRare results for this sample.")
      }
    }

    if(currentSample %in% resObj$whichSamplesMayHaveRareSigs & !forceCommon){
      highestCosSimSig <- resObj$candidateRareSigs[[currentSample]][which.max(resObj$candidateRareSigsCosSim[[currentSample]])]
      if(!is.null(forceRareSig)) highestCosSimSig <- forceRareSig
      rareSigChoice[[currentSample]] <- highestCosSimSig
      selectedExp <- t(resObj$samples[[currentSample]]$fitWithRare[[highestCosSimSig]]$exposures)
      exposures_merge[currentSample,rownames(selectedExp)] <- selectedExp
      # exposures_merge[currentSample,"unassigned"] <- resObj$samples[[currentSample]]$fitWithRare[[highestCosSimSig]]$unassigned_muts
      # collect bootstraps
      if(resObj$useBootstrap) bootstrap_exposures_samples[[currentSample]] <- resObj$samples[[currentSample]]$fitWithRare[[highestCosSimSig]]$bootstrap_exposures_samples[[1]]
    }else{
      if(forceCommon) rareSigChoice[[currentSample]] <- "commonOnly"
      selectedExp <- t(resObj$samples[[currentSample]]$fitCommonOnly$exposures)
      exposures_merge[currentSample,rownames(selectedExp)] <- selectedExp
      # exposures_merge[currentSample,"unassigned"] <- resObj$samples[[currentSample]]$fitCommonOnly$unassigned_muts
      # collect bootstraps
      if(resObj$useBootstrap) bootstrap_exposures_samples[[currentSample]] <- resObj$samples[[currentSample]]$fitCommonOnly$bootstrap_exposures_samples[[1]]
    }
  }

  resObj$exposures <- exposures_merge[,apply(exposures_merge, 2, sum)>0,drop=F]
  resObj$rareSigChoice <- rareSigChoice
  if(resObj$useBootstrap){
    # save the bootstraps organised by samples
    resObj$bootstrap_exposures_samples <- bootstrap_exposures_samples
    # also save the bootstraps organised by bootstrap
    nboots <- ncol(bootstrap_exposures_samples[[1]])
    sample_names <- names(bootstrap_exposures_samples)
    sigslist <- c()
    for(sn in sample_names) sigslist <- union(sigslist,rownames(bootstrap_exposures_samples[[sn]]))
    bootstrap_exposures <- list()
    for(i in 1:nboots){
      current_subs <- matrix(0,nrow = length(sigslist),ncol = length(sample_names),dimnames = list(sigslist,sample_names))
      for(sn in sample_names) {
        current_subs[rownames(bootstrap_exposures_samples[[sn]]),sn] <- bootstrap_exposures_samples[[sn]][,i]
      }
      bootstrap_exposures[[i]] <- current_subs
    }
    resObj$bootstrap_exposures <- bootstrap_exposures
  }else{
    resObj$bootstrap_exposures_samples <- NULL
  }

  return(resObj)
}


flexconstr_sigfit <- function(P,m, mut_tol=0, allmut_tolratio=0.01){
  G <- -P
  H <- -m
  # test to flexibilise the constraints on the residual part
  allmut <-sum(m)
  H<- -m-mut_tol*m-allmut_tolratio*allmut

  # G<-NULL
  # H<-NULL

  # add positivity requirement for exposures e >=0
  if (ncol(P)>1){
    coef_nn_e<-diag(1,ncol(P))
    rside_nn_e <- rep(0,ncol(P))
    # summarising both constraints into matrices G and H
    G <- rbind(G,coef_nn_e)
    H <- c(H,rside_nn_e)
  }else{
    G <- c(G,1)
    H <- c(H,0)
  }

  # solving least squares on (Pe-m)^2  with the constraints (-Pe>=-m-mut_tol-allmut_tolratio*) and (e >=0)  => Gx>=h
  e <- limSolve::lsei(A = P, B= m, G = G,
                      H = H, E=NULL, F=NULL, Wx = NULL, Wa = NULL, type = 2, tol = sqrt(.Machine$double.eps),
                      tolrank = NULL, fulloutput = TRUE, verbose = TRUE)

  return(e)

}

flexconstr_sigfit_multipleSamples <- function(P,M, mut_tol=0, allmut_tolratio=0.01){
  E <- list()
  for(i in 1:ncol(M)){
    m <- M[,i,drop=F]
    e <- flexconstr_sigfit(P,m, mut_tol, allmut_tolratio)
    E[[colnames(M)[i]]] <- as.matrix(e$X)
  }
  E <- do.call(cbind,E)
  colnames(E) <- colnames(M)
  return(E)
}

partialNMFresidual <- function(catalogues,commonSignatures,E,max.iter = 100){
  Rlist <- list()
  for (i in 1:ncol(E)){
    # i <- 1
    res.nnmf <- NNLM::nnmf(as.matrix(catalogues[,i,drop=F]),
                           init = list(W0 = as.matrix(commonSignatures),H1 = as.matrix(E[,i,drop=F])),
                           loss = "mkl",
                           method = "lee",
                           k = 1,
                           max.iter = max.iter,check.k = FALSE,show.warning = F,verbose = F)
    currentR <- apply(res.nnmf$W[,1,drop=F],2,function(x)x/sum(x))
    colnames(currentR) <- colnames(catalogues)[i]
    Rlist[[i]] <- currentR
  }
  return(do.call(cbind,Rlist))
}

#' @export
giniCoeff <- function(x){
  meanx <- mean(x)
  sumxdiff <- 0
  for (i in x) {
    for (j in x){
      sumxdiff <- sumxdiff + abs(i-j)
    }
  }
  return((sumxdiff)/(2*(length(x)^2)*meanx))
}



MAD <- function(a,b){
  mean(abs(as.matrix(a)-as.matrix(b)))
}

NMAD <- function(a,b){
  mean(abs(as.matrix(a)-as.matrix(b)))/sum(as.matrix(a))
}

# dataMatrix can be either a catalogue or signature matrix
getTypeOfMutationsFromChannels <- function(dataMatrix){
  SBSchannels <- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T",
                   "A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T",
                   "A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T",
                   "A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T",
                   "A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T",
                   "A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")

  DBSchannelsZou <- c("AA>CC","AA>CG","AA>CT","AA>GC","AA>GG","AA>GT","AA>TC","AA>TG","AA>TT",
                        "AC>CA","AC>CG","AC>CT","AC>GA","AC>GG","AC>GT","AC>TA","AC>TG","AC>TT",
                        "AG>CA","AG>CC","AG>CT","AG>GA","AG>GC","AG>GT","AG>TA","AG>TC","AG>TT",
                        "AT>CA","AT>CC","AT>CG","AT>GA","AT>GC","AT>TA",
                        "CA>AC","CA>AG","CA>AT","CA>GC","CA>GG","CA>GT","CA>TC","CA>TG","CA>TT",
                        "CC>AA","CC>AG","CC>AT","CC>GA","CC>GG","CC>GT","CC>TA","CC>TG","CC>TT",
                        "CG>AA","CG>AC","CG>AT","CG>GA","CG>GC","CG>TA",
                        "GA>AC","GA>AG","GA>AT","GA>CC","GA>CG","GA>CT","GA>TC","GA>TG","GA>TT",
                        "GC>AA","GC>AG","GC>AT","GC>CA","GC>CG","GC>TA",
                        "TA>AC","TA>AG","TA>AT","TA>CC","TA>CG","TA>GC")
  DBSchannelsAlexandrov <-    c("AC>CA","AC>CG","AC>CT","AC>GA","AC>GG","AC>GT","AC>TA","AC>TG","AC>TT",
                                  "AT>CA","AT>CC","AT>CG","AT>GA","AT>GC","AT>TA",
                                  "CC>AA","CC>AG","CC>AT","CC>GA","CC>GG","CC>GT","CC>TA","CC>TG","CC>TT",
                                  "CG>AT","CG>GC","CG>GT","CG>TA","CG>TC","CG>TT",
                                  "CT>AA","CT>AC","CT>AG","CT>GA","CT>GC","CT>GG","CT>TA","CT>TC","CT>TG",
                                  "GC>AA","GC>AG","GC>AT","GC>CA","GC>CG","GC>TA",
                                  "TA>AT","TA>CG","TA>CT","TA>GC","TA>GG","TA>GT",
                                  "TC>AA","TC>AG","TC>AT","TC>CA","TC>CG","TC>CT","TC>GA","TC>GG","TC>GT",
                                  "TG>AA","TG>AC","TG>AT","TG>CA","TG>CC","TG>CT","TG>GA","TG>GC","TG>GT",
                                  "TT>AA","TT>AC","TT>AG","TT>CA","TT>CC","TT>CG","TT>GA","TT>GC","TT>GG")
  SVchannels <- c('clustered_del_1-10Kb', 'clustered_del_10-100Kb', 'clustered_del_100Kb-1Mb', 'clustered_del_1Mb-10Mb', 'clustered_del_>10Mb',
                  'clustered_tds_1-10Kb', 'clustered_tds_10-100Kb', 'clustered_tds_100Kb-1Mb', 'clustered_tds_1Mb-10Mb', 'clustered_tds_>10Mb',
                  'clustered_inv_1-10Kb', 'clustered_inv_10-100Kb', 'clustered_inv_100Kb-1Mb', 'clustered_inv_1Mb-10Mb', 'clustered_inv_>10Mb',
                  'clustered_trans',
                  'non-clustered_del_1-10Kb', 'non-clustered_del_10-100Kb', 'non-clustered_del_100Kb-1Mb', 'non-clustered_del_1Mb-10Mb', 'non-clustered_del_>10Mb',
                  'non-clustered_tds_1-10Kb', 'non-clustered_tds_10-100Kb', 'non-clustered_tds_100Kb-1Mb', 'non-clustered_tds_1Mb-10Mb', 'non-clustered_tds_>10Mb',
                  'non-clustered_inv_1-10Kb', 'non-clustered_inv_10-100Kb', 'non-clustered_inv_100Kb-1Mb', 'non-clustered_inv_1Mb-10Mb', 'non-clustered_inv_>10Mb',
                  'non-clustered_trans')

  muttype <- "generic"

  if(nrow(dataMatrix)==length(SBSchannels)){
    if(all(rownames(dataMatrix)==SBSchannels)) muttype <- "subs"
  }else if(nrow(dataMatrix)==length(DBSchannelsZou)){
    if(all(rownames(dataMatrix)==DBSchannelsZou)) {
      muttype <- "DNV"
    }else if(all(rownames(dataMatrix)==DBSchannelsAlexandrov)){
      muttype <- "DNV"
    }
  }else if(nrow(dataMatrix)==length(SVchannels)){
    if(all(rownames(dataMatrix)==SVchannels)) muttype <- "rearr"
  }

  return(muttype)
}

# plotSignatures wrapper function, it will determine the type of signatures and plot using the appropriate function
#' Plot Signatures with automated detection of type of mutations
#'
#' This function checks the channels of the input matrix, determines the type of mutations and plots using the
#' most appropriate signature plot.
#'
#' @param signature_data_matrix matrix of signatures, signatures as columns and channels as rows
#' @param output_file set output file, should end with ".jpg" or ".pdf". If output_file==null, output will not be to a file, but will still run the plot functions. The option output_file==null can be used to add this plot to a larger output file.
#' @param plot_sum whether the sum of the channels should be plotted. If plotting signatures this should be FALSE, but if plotting sample catalogues, this can be set to TRUE to display the number of mutations in each sample.
#' @param overall_title set the overall title of the plot
#' @param mar set the option par(mar=mar)
#' @param howManyInOnePage how many signatures or catalogues should be plotted on one page. Multiple pages are plotted if more signatures/catalogues to plot have been requested
#' @param ncolumns how many columns should be used to arrange the signatures/catalogues to plot
#' @export
plotSignatures <- function(signature_data_matrix,
                           output_file = NULL,
                           plot_sum = TRUE,
                           overall_title = "",
                           add_to_titles = NULL,
                           mar=NULL,
                           howManyInOnePage=100,
                           ncolumns=1){
  # identify the type of mutations
  typeofmuts <- getTypeOfMutationsFromChannels(signature_data_matrix)
  if(typeofmuts=="subs"){
    plotSubsSignatures(signature_data_matrix = signature_data_matrix,
                       output_file = output_file,
                       plot_sum = plot_sum,
                       overall_title = overall_title,
                       add_to_titles = add_to_titles,
                       mar = mar,
                       howManyInOnePage = howManyInOnePage,
                       ncolumns = ncolumns)
  }else if(typeofmuts=="rearr"){
    plotRearrSignatures(signature_data_matrix = signature_data_matrix,
                        output_file = output_file,
                        plot_sum = plot_sum,
                        overall_title = overall_title,
                        add_to_titles = add_to_titles,
                        mar = mar,
                        howManyInOnePage = howManyInOnePage,
                        ncolumns = ncolumns)
  }else if(typeofmuts=="DNV"){
    plotDNVSignatures(signature_data_matrix = signature_data_matrix,
                      output_file = output_file,
                      plot_sum = plot_sum,
                      overall_title = overall_title,
                      add_to_titles = add_to_titles,
                      mar = mar,
                      howManyInOnePage = howManyInOnePage,
                      ncolumns = ncolumns)
  }else{
    plotGenericSignatures(signature_data_matrix = signature_data_matrix,
                          output_file = output_file,
                          plot_sum = plot_sum,
                          overall_title = overall_title,
                          add_to_titles = add_to_titles,
                          mar = mar,
                          howManyInOnePage = howManyInOnePage,
                          ncolumns = ncolumns)
  }
}

drawCircle <- function(radius=1,
                       position=c(0,0),
                       col=NA,         #this is a polygon parameter
                       border=NULL){   #this is a polygon parameter

  displCoor <- position


  arcLines <- getArcCoordinates(0,2*pi,radius = radius,nSegmentsForEachRadiant = 20)
  polygon(displCoor[1]+c(arcLines$x),
          displCoor[2]+c(arcLines$y),
          col = col,
          border = border)
}

polarToCartesian <- function(angle,
                             radius=1){
  c(radius*sin(angle),radius*cos(angle))
}

getArcCoordinates <- function(fromAngle,
                              toAngle,
                              radius=1,
                              nSegmentsForEachRadiant=100){
  if(fromAngle>toAngle) toAngle <- toAngle + 2*pi
  angleToTravel <- toAngle - fromAngle
  nsegments <- ceiling(angleToTravel*nSegmentsForEachRadiant)
  angles <- seq(fromAngle,toAngle,length.out = nsegments)
  resPos <- sapply(angles, function(x){
    polarToCartesian(x,radius)
  })
  resList <- list()
  resList$x <- resPos[1,]
  resList$y <- resPos[2,]
  return(resList)
}

#' Plot Matrix
#'
#' This function plots a matrix of values. Data is visualised as circles scaled with respect to the
#' largest value in the matrix, with the actual number shown. If thresholdMark is specified, then a
#' different colour is used for entries that are equal or above the threshold.
#'
#' @param dataMatrix a data matrix or data frame
#' @param output_file if an output file name is given (must be pdf), the matrix will be plotted to file
#' @param thresholdMark threshold for using a different colour to highlight entries above a threshold
#' @param ndigitsafterzero specify how many digits after the zero should be used to show the actual numbers
#' @param cex.numbers scale the text used for the numbers in the matrix
#' @param circlesColBasic colour used for the circles
#' @param circlesColHighlight colour used for the circles that pass the thresholdMark
#' @export
plotMatrix <- function(dataMatrix,
                       output_file = NULL,
                       thresholdMark = NULL,
                       ndigitsafterzero = 2,
                       cex.numbers = 0.7,
                       circlesColBasic = "#A1CAF1",
                       circlesColHighlight = "#F6A600"){

  maxncharSigs <- max(sapply(rownames(dataMatrix),nchar))
  maxncharSamples <- max(sapply(colnames(dataMatrix),nchar))
  mar1 <- 0.6*maxncharSigs+1.2
  mar2 <- 0.6*maxncharSamples+1.2
  mar3 <- 2
  mar4 <- 2
  width <- 0.5 + 0.3*nrow(dataMatrix) + 0.125*maxncharSamples
  height <- 0.5 + 0.3*ncol(dataMatrix) + 0.125*maxncharSigs
  if(!is.null(output_file)) cairo_pdf(filename = output_file,width = width,height = height)
  par(mfrow=c(1,1))
  par(mar=c(mar1,mar2,mar3,mar4))

  plot(1, type="n", xlab="", ylab="", xlim=c(0.5,nrow(dataMatrix)+0.5), ylim=c(ncol(dataMatrix)+0.5,0.5),
       xaxt = 'n', yaxt = 'n',bty = 'n',xaxs="i",yaxs="i")
  abline(h=1:ncol(dataMatrix),lty=3,col="lightgrey",lwd=3)
  abline(v=1:nrow(dataMatrix),lty=3,col="lightgrey",lwd=3)
  axis(1,labels = rownames(dataMatrix),at = 1:nrow(dataMatrix),las=2,lwd = 0,lwd.ticks = 1,cex.axis = 1)
  axis(2,labels = colnames(dataMatrix),at = 1:ncol(dataMatrix),las=2,lwd = 0,lwd.ticks = 1,cex.axis = 1)

  toPlot <- dataMatrix
  for(i in 1:ncol(dataMatrix)) toPlot[,i] <- sprintf(paste0("%.",ndigitsafterzero,"f"),dataMatrix[,i])
  toPlot[toPlot=="0" | toPlot=="-0"] <- ""

  circleDim <- dataMatrix/max(dataMatrix)*5

  for(i in 1:ncol(circleDim)) {
    for(j in 1:nrow(circleDim)){
      usecol <- circlesColBasic
      if(!is.null(thresholdMark)) {
        if(thresholdMark <= dataMatrix[j,i]) usecol <- circlesColHighlight
      }
      drawCircle(radius = circleDim[j,i]/10,position = c(j,i),col = usecol,border = NA)
    }
  }
  for(i in 1:ncol(toPlot)) text(y = rep(i,nrow(toPlot)), x = 1:nrow(toPlot),labels = toPlot[,i],cex = cex.numbers)
  if(!is.null(output_file)) dev.off()
}

#' Plot the results from the Fit or FitMS function
#'
#' Plotting of the results obtained with the Fit or FitMS function. Output adapts based on the options used in during fitting.
#' You can use this plot function instead of using plotFit or plotFitMS. The function plotFitResults will infer whether Fit or FitMS
#' was used and use the appropriate plot function.
#'
#' @param fitObj object obtained from the Fit or FitMS function
#' @param outdir output directory where the results should be saved/plotted
#' @export
#' @examples
#' res <- Fit(catalogues,getOrganSignatures("Breast"))
#' plotFitResults(res,"results/")
plotFitResults <- function(fitObj,
                           outdir = ""){
  if(!is.null(fitObj$fitAlgorithm)){
    if(fitObj$fitAlgorithm=="Fit"){
      plotFit(fitObj = fitObj,
              outdir = outdir)
    }else if(fitObj$fitAlgorithm=="FitMS"){
      plotFitMS(fitMSobj = fitObj,
                outdir = outdir)
    }else{
      message("[error plotFitResults] unknown fitAlgorithm attribute, expecting Fit or FitMS. Input fitObj not recognised. ")
    }
  }else{
    message("[error plotFitResults] missing fitAlgorithm attribute, fitObj not recognised.")
  }
}

#' Plot the results from the Fit function
#'
#' Plotting of the results obtained with the Fit function. Output adapts based on the options used in the Fit function
#'
#' @param fitObj object obtained from the Fit function
#' @param outdir output directory where the results should be saved/plotted
#' @param samplesInSubdir if TRUE, move the sample files to a subdirectory named "samples"
#' @export
#' @examples
#' res <- Fit(catalogues,getOrganSignatures("Breast"))
#' plotFit(res,"results/")
plotFit <- function(fitObj,
                      outdir = "",
                      samplesInSubdir = TRUE){

  # some checks on outdir
  if(is.null(outdir)) {
    message("[error plotFit] please specify outdir")
    return(NULL)
  }
  if(outdir != ""){
    if(substr(outdir,nchar(outdir),nchar(outdir)) != "/") outdir <- paste0(outdir,"/")
  }

  # now let's plot
  if(outdir != "") dir.create(outdir,showWarnings = F,recursive = T)

  #function to draw a legend for the heatmap of the correlation matrix
  draw_legend <- function(col,xl,xr,yb,yt,textx){
    par(xpd=TRUE)
    rect(xl,yb,xr,yt)
    rect(
      xl,
      head(seq(yb,yt,(yt-yb)/length(col)),-1),
      xr,
      tail(seq(yb,yt,(yt-yb)/length(col)),-1),
      col=col,border = NA
    )
    text(x = textx, y = yt,labels = "1")
    text(x = textx, y = (yt-yb)/2,labels = "0")
    text(x = textx, y = yb,labels = "-1")
    par(xpd=FALSE)
  }

  reconstructed <- round(as.matrix(fitObj$signatures) %*% t(fitObj$exposures[,1:(ncol(fitObj$exposures)-1),drop=F]))

  #plot and save exposures
  sums_exp <- apply(fitObj$exposures,1,sum)
  denominator <- matrix(sums_exp,nrow = nrow(fitObj$exposures),ncol = ncol(fitObj$exposures),byrow = FALSE)
  exposuresProp <- (fitObj$exposures/denominator*100)
  # case of empty catalogues
  exposuresProp[sums_exp==0,] <- 0

  # #plot and save exposures
  # sums_exp <- apply(fitObj$catalogues, 2, sum)
  # exposures <- rbind(fitObj$exposures,fitObj$unassigned_muts)
  # rownames(exposures)[nrow(exposures)] <- "unassigned"
  # denominator <- matrix(sums_exp,nrow = nrow(exposures),ncol = ncol(exposures),byrow = TRUE)
  # exposuresProp <- (exposures/denominator*100)
  # # case of empty catalogues
  # exposuresProp[,sums_exp==0] <- 0

  file_table_exp <- paste0(outdir,"exposures.tsv")
  # change to pdf later
  file_plot_exp <- paste0(outdir,"exposures.pdf")
  file_plot_expProp <- paste0(outdir,"exposures_prop.pdf")
  plotMatrix(as.data.frame(t(fitObj$exposures)),output_file = file_plot_exp,ndigitsafterzero = 0)
  plotMatrix(as.data.frame(t(exposuresProp)),output_file = file_plot_expProp,ndigitsafterzero = 0)

  write.table(fitObj$exposures,file = file_table_exp,
              sep = "\t",col.names = TRUE,row.names = TRUE,quote = FALSE)

  #provide a series of plots for each sample
  if(samplesInSubdir){
    outdir <- paste0(outdir,"samples/")
    dir.create(outdir,recursive = TRUE,showWarnings = FALSE)
  }
  howmanyplots <- ncol(fitObj$catalogues)
  for(p in 1:howmanyplots){
    # p <- 1
    currentSample <- colnames(fitObj$catalogues)[p]
    fitIsEmpty <- sum(fitObj$exposures[p,])==0
    unassigned_mut <- sprintf("%.2f",(fitObj$unassigned_muts_perc[p]))
    cos_sim <- sprintf("%.2f",cos_sim(fitObj$catalogues[,p,drop=FALSE],reconstructed[,p,drop=FALSE]))

    sigMatrix <- fitObj$catalogues[,p,drop=FALSE]
    addToTitle <- c("Catalogue")
    if(!fitIsEmpty){
      sigMatrix <- cbind(sigMatrix,reconstructed[,p,drop=FALSE])
      sigMatrix <- cbind(sigMatrix,fitObj$catalogues[,p,drop=FALSE] - reconstructed[,p,drop=FALSE])
      addToTitle <- c(addToTitle,
                      paste0("Model (CosSim ",cos_sim,")"),
                      paste0("Difference (Unassigned ",unassigned_mut,"%)"))
    }
    colnames(sigMatrix) <- paste0(colnames(sigMatrix)," ",addToTitle)

    plotSignatures(signature_data_matrix = sigMatrix,add_to_title = NULL,
                   output_file = paste0(outdir,"signatureFit_",p,"of",howmanyplots,"_",currentSample,"_pointEstimate.pdf"))

    if(fitObj$useBootstrap & !fitIsEmpty){
      #4 bootstraps
      consensusCol <- "#BE0032"
      thresholdCol <- "#8DB600"
      bootsCol <- "#A1CAF1"

      thresholdText <- ""
      if(fitObj$exposureFilterType=="fixedThreshold"){
        thresholdText <- paste0("threshold=",fitObj$threshold_percent,"%")
      }else if(fitObj$exposureFilterType=="giniScaledThreshold") {
        thresholdText <- paste0("threshold=(1-Gini)*",fitObj$giniThresholdScaling)
      }

      maxnchar <- max(sapply(colnames(fitObj$signatures),nchar))
      mar1 <- 0.6*maxnchar+1.2
      mar2 <- 4
      mar3 <- 5
      mar4 <- 2
      width <- max(4.5,0.5*ncol(fitObj$signatures)+1)
      height <- 3 + 0.125*maxnchar
      cairo_pdf(filename = paste0(outdir,"signatureFit_",p,"of",howmanyplots,"_",currentSample,"_Bootstrap.pdf"),width = width,height = height)
      par(mfrow=c(1,1))
      par(mar=c(mar1,mar2,mar3,mar4))
      boxplot(t(fitObj$bootstrap_exposures_samples[[p]]),las=3,cex.axes=0.9,
              ylab="n mutations",
              ylim=c(0,max(fitObj$bootstrap_exposures_samples[[p]])),cex.main = 0.9,border="#848482",
              main=paste0("Exposures, of ",rownames(fitObj$exposures)[p],"\n",thresholdText,", p-value=",fitObj$threshold_p.value,", n=",fitObj$nboot))
      if(fitObj$exposureFilterType=="fixedThreshold"){
        abline(h=fitObj$threshold_percent/100*sum(fitObj$catalogues[,p,drop=FALSE]),col=thresholdCol,lwd = 2)
      }else if(fitObj$exposureFilterType=="giniScaledThreshold") {
        sigInvGini <- 1 - apply(fitObj$signatures,2,giniCoeff)
        giniThresholdPerc <- fitObj$giniThresholdScaling*sigInvGini
        giniThreshold <- giniThresholdPerc/100*sum(fitObj$catalogues[,p,drop=FALSE])
        for(si in 1:length(giniThreshold)) lines(x = c(si-0.5,si+0.5),y = rep(giniThreshold[si],2),col=thresholdCol,lwd = 2)
      }
      points(1:(length(fitObj$exposures[p,])-1),fitObj$exposures[p,1:(ncol(fitObj$exposures)-1)],col=consensusCol,pch = 16)
      legend(x="topleft",legend = c("consensus exposures"),col = consensusCol,pch = 16,cex = 0.9,bty = "n",inset = c(0,-0.14),xpd = TRUE)
      legend(x="topright",legend = c("threshold"),col = thresholdCol,lty = 1,cex = 0.9,bty = "n",inset = c(0,-0.14),xpd = TRUE,lwd = 2)
      dev.off()

      if(ncol(fitObj$signatures)>1){
        #5 top correlated signatures
        res.cor <- suppressWarnings(cor(t(fitObj$bootstrap_exposures_samples[[p]]),method = "spearman"))
        res.cor_triangular <- res.cor
        res.cor_triangular[row(res.cor)+(ncol(res.cor)-col(res.cor))>=ncol(res.cor)] <- 0
        res.cor_triangular_label <- matrix(sprintf("%0.2f",res.cor_triangular),nrow = nrow(res.cor_triangular))
        res.cor_triangular_label[row(res.cor)+(ncol(res.cor)-col(res.cor))>=ncol(res.cor)] <- ""

        mar1 <- 0.6*maxnchar+1.2
        mar2 <- 0.6*maxnchar+1.2
        mar3 <- 4
        mar4 <- 4
        width <- 1 + 0.3*ncol(fitObj$signatures) + 0.125*maxnchar
        height <- 1 + 0.3*ncol(fitObj$signatures) + 0.125*maxnchar
        cairo_pdf(filename = paste0(outdir,"signatureFit_",p,"of",howmanyplots,"_",currentSample,"_Bootstrap_Correlation.pdf"),width = width,height = height)
        par(mfrow=c(1,1))
        par(mar=c(mar1,mar2,mar3,mar4))
        par(xpd=FALSE)
        col<- colorRampPalette(c("blue", "white", "red"))(51)
        image(res.cor_triangular,col = col,zlim = c(-1,1), axes=F,main="Exposures Correlation\n(spearman)")
        extrabit <- 1/(ncol(fitObj$signatures)-1)/2
        abline(h=seq(0-extrabit,1+extrabit,length.out = ncol(fitObj$signatures)+1),col="grey",lty=2)
        abline(v=seq(0-extrabit,1+extrabit,length.out = ncol(fitObj$signatures)+1),col="grey",lty=2)
        axis(2,at = seq(0,1,length.out = ncol(fitObj$signatures)),labels = colnames(fitObj$signatures),las=1,cex.lab=0.8)
        axis(1,at = seq(0,1,length.out = ncol(fitObj$signatures)),labels = colnames(fitObj$signatures),las=2,cex.lab=0.8)
        # draw_legend(col,1.25,1.3,0,1,1.2)
        draw_legend(col,1+4*extrabit,1+4.5*extrabit,0,1,1+3*extrabit)
        dev.off()

        #6 some correlation plots
        howmanycorrtoplot <- min(3,ncol(fitObj$signatures)*(ncol(fitObj$signatures)-1)/2)
        cairo_pdf(filename = paste0(outdir,"signatureFit_",p,"of",howmanyplots,"_",currentSample,"_Bootstrap_CorrelationExamples.pdf"),width = 2.2*howmanycorrtoplot,height = 2.5)
        par(mfrow=c(1,howmanycorrtoplot))
        vals <- res.cor_triangular[order(abs(res.cor_triangular),decreasing = TRUE)]
        for (j in 1:howmanycorrtoplot){
          pos <- which(vals[j]==res.cor_triangular,arr.ind = TRUE)
          mainpar <- paste0("Exposures across bootstraps, n=",fitObj$nboot,"\nspearman correlation ",sprintf("%.2f",vals[j]))
          plot(fitObj$bootstrap_exposures_samples[[p]][pos[1],],fitObj$bootstrap_exposures_samples[[p]][pos[2],],
               xlab = colnames(fitObj$signatures)[pos[1]],
               ylab = colnames(fitObj$signatures)[pos[2]],
               main=mainpar,col=bootsCol,pch = 16,cex.main = 0.9)

        }
        dev.off()
      }

    }
  }

}

#' Plot the results from the FitMS function
#'
#' Plotting of the results obtained with the FitMS function. Output adapts based on the options used in the FitMS function
#'
#' @param fitObj object obtained from the FitMS function
#' @param outdir output directory where the results should be saved/plotted
#' @export
#' @examples
#' res <- FitMS(catalogues,"Breast")
#' plotFitMS(res,"results/")
plotFitMS <- function(fitMSobj,
                      outdir = ""){
  # some checks on outdir
  if(is.null(outdir)) {
    message("[error plotFitMS] please specify outdir")
    return(NULL)
  }
  if(outdir != ""){
    if(substr(outdir,nchar(outdir),nchar(outdir)) != "/") outdir <- paste0(outdir,"/")
  }
  # now let's plot
  if(outdir != "") dir.create(outdir,showWarnings = F,recursive = T)

  # now plot some generic info
  if(length(fitMSobj$whichSamplesMayHaveRareSigs)>0){
    summaryFits <- NULL
    rowi <- 1
    for (i in 1:length(fitMSobj$whichSamplesMayHaveRareSigs)){
      # i <- 1
      currentSample <- fitMSobj$whichSamplesMayHaveRareSigs[i]
      for (j in 1:length(fitMSobj$candidateRareSigs[[currentSample]])){
        summaryFits <- rbind(summaryFits,data.frame(sample=currentSample,
                                                    candidate=fitMSobj$candidateRareSigs[[currentSample]][j],
                                                    cossim=fitMSobj$candidateRareSigsCosSim[[currentSample]][j],
                                                    selected=fitMSobj$rareSigChoice[[currentSample]]==fitMSobj$candidateRareSigs[[currentSample]][j],
                                                    stringsAsFactors = F,row.names = rowi))
        rowi <- rowi+1
      }
    }
    write.table(summaryFits,file = paste0(outdir,"candidateRareSigs.tsv"),
                sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
  }
  #plot and save exposures
  sums_exp <- apply(fitMSobj$exposures,1,sum)
  denominator <- matrix(sums_exp,nrow = nrow(fitMSobj$exposures),ncol = ncol(fitMSobj$exposures),byrow = FALSE)
  exposuresProp <- (fitMSobj$exposures/denominator*100)
  # case of empty catalogues
  exposuresProp[sums_exp==0,] <- 0

  file_plot_exp <- paste0(outdir,"exposures.pdf")
  file_plot_expProp <- paste0(outdir,"exposures_prop.pdf")
  file_table_exp <- paste0(outdir,"exposures.tsv")
  plotMatrix(as.data.frame(t(fitMSobj$exposures)),output_file = file_plot_exp,ndigitsafterzero = 0)
  plotMatrix(as.data.frame(t(exposuresProp)),output_file = file_plot_expProp,ndigitsafterzero = 0)
  write.table(fitMSobj$exposures,file = file_table_exp,
              sep = "\t",col.names = TRUE,row.names = TRUE,quote = FALSE)

  # now for each sample plot the chosen solution and the alternative solutions
  for (s in colnames(fitMSobj$catalogues)){
    # s <- colnames(fitMSobj$catalogues)[1]
    if(s %in% fitMSobj$whichSamplesMayHaveRareSigs){
      currentOutDir <- paste0(outdir,"selectedSolution_",s,"/")
      selectedSolution <- fitMSobj$rareSigChoice[[s]]
      otherSolutions <- setdiff(fitMSobj$candidateRareSigs[[s]],selectedSolution)
      if(selectedSolution=="commonOnly"){
        plotFit(fitMSobj$samples[[s]]$fitCommonOnly,outdir = currentOutDir,samplesInSubdir = F)
      }else{
        plotFit(fitMSobj$samples[[s]]$fitWithRare[[selectedSolution]],outdir = currentOutDir,samplesInSubdir = F)
      }

      # plot other candidates
      currentOutDir <- paste0(outdir,"otherSolutions_",s,"/commonOnly/")
      if(!selectedSolution=="commonOnly") plotFit(fitMSobj$samples[[s]]$fitCommonOnly,outdir = currentOutDir,samplesInSubdir = F)

      if(length(otherSolutions)>0){
        for (i in 1:length(otherSolutions)){
          # i <- 1
          currentOutDir <- paste0(outdir,"otherSolutions_",s,"/rareSigs_",i,"/")
          plotFit(fitMSobj$samples[[s]]$fitWithRare[[otherSolutions[i]]],outdir = currentOutDir,samplesInSubdir = F)
        }
      }

    }else{
      currentOutDir <- paste0(outdir,"selectedSolution_",s,"/")
      plotFit(fitMSobj$samples[[s]]$fitCommonOnly,outdir = currentOutDir,samplesInSubdir = F)
    }
  }
}

vectorToJSON <- function(v,
                         vectorname,
                         nindent=0){
  indent <- rep("\t",nindent)
  cat(indent,"\"",vectorname,"\": ",sep = "")
  if(!is.null(names(v))){
    cat(indent,"{\n",sep = "")
    for (i in 1:length(v)) {
      cat(indent,"\t\"",names(v)[i],"\": ",sep = "")
      if(typeof(v[i])=="double"){
        cat(v[i],sep = "")
      }else if (typeof(v[i])=="character"){
        cat("\"",v[i],"\"",sep = "")
      }
      if(i<length(v)) cat(",")
      cat("\n")
    }
    cat(indent,"}",sep = "")
  }else{
    cat(indent,"[",sep = "")
    for (i in 1:length(v)) {
      if(typeof(v[i])=="double"){
        cat(v[i],sep = "")
      }else if (typeof(v[i])=="character"){
        cat("\"",v[i],"\"",sep = "")
      }
      if(i<length(v)) cat(", ")
    }
    cat("]",sep = "")
  }
}

#internal fuction for plotting a table in JSON
tableToJSON <- function(tab,
                        tablename,
                        nindent = 0){
  indent <- rep("\t",nindent)
  cat(indent,"\"",tablename,"\": {\n",sep = "")
  for(i in 1:ncol(tab)){
    sname <- colnames(tab)[i]
    cat(indent,"\t\"",sname,"\": {\n",sep = "")

    for(j in 1:nrow(tab)){
      rname <- rownames(tab)[j]
      cat(indent,"\t\t\"",rname,"\": ",tab[j,i],sep = "")
      if(j<nrow(tab)){
        cat(",\n")
      }else{
        cat("\n")
      }
    }

    cat(indent,"\t}",sep = "")
    if(i<ncol(tab)){
      cat(",\n")
    }else{
      cat("\n")
    }
  }
  cat(indent,"}",sep = "")
}


#' Export the results from the Fit function to JSON
#'
#' Exports of the results obtained with the Fit function to JSON. Output adapts based on the options used in the Fit function
#'
#' @param fitObj object obtained from the Fit function
#' @param filename file name where to save the JSON string
#' @param nindent number of tabs of indentation to be used, default=0. This is useful in case one wants to embed the output in a larger JSON document
#' @param disableFileWrite if TRUE the JSON lines are simply printed to screen. This option is used when fitToJSON is called from another function like fitMStoJSON, which itself opens a sink() to write to.
#' @param blockname if specified, the first line will show the given name as the block name, which is useful if the output is embedded in a larger JSON file
#' @export
#' @examples
#' res <- Fit(catalogues,getOrganSignatures("Breast"))
#' fitToJSON(res,"export.json")
fitToJSON <- function(fitObj,
                      filename = "export.json",
                      nindent = 0,
                      disableFileWrite = FALSE,
                      blockname = NULL){
  #plot consensus exposures file with metadata
  indent <- rep("\t",nindent)
  if(!disableFileWrite) sink(file = filename)
  if(is.null(blockname)){
    cat(indent,"{\n",sep = "")
  }else{
    cat(indent,"\"",blockname,"\": {\n",sep = "")
  }
  cat(indent,"\t\"fitAlgorithm\": ",fitObj$fitAlgorithm,",\n",sep = "")
  cat(indent,"\t\"nsamples\": ",nrow(fitObj$exposures),",\n",sep = "")
  cat(indent,"\t\"nsignatures\": ",ncol(fitObj$signatures),",\n",sep = "")
  cat(indent,"\t\"nchannels\": ",nrow(fitObj$catalogues),",\n",sep = "")
  cat(indent,"\t\"method\": \"",fitObj$method,"\",\n",sep = "")
  cat(indent,"\t\"exposureFilterType\": \"",fitObj$exposureFilterType,"\",\n",sep = "")
  if(!is.na(fitObj$threshold_percent)){
    cat(indent,"\t\"threshold_percent\": ",fitObj$threshold_percent,",\n",sep = "")
  }else{
    cat(indent,"\t\"threshold_percent\": null,\n",sep = "")
  }
  if(!is.na(fitObj$giniThresholdScaling)){
    cat(indent,"\t\"giniThresholdScaling\": ",fitObj$giniThresholdScaling,",\n",sep = "")
  }else{
    cat(indent,"\t\"giniThresholdScaling\": null,\n",sep = "")
  }
  cat(indent,"\t\"useBootstrap\": ",ifelse(fitObj$useBootstrap,"true","false"),",\n",sep = "")
  if(!is.na(fitObj$nboot)){
    cat(indent,"\t\"nboot\": ",fitObj$nboot,",\n",sep = "")
  }else{
    cat(indent,"\t\"nboot\": null,\n",sep = "")
  }
  if(!is.na(fitObj$threshold_p.value)){
    cat(indent,"\t\"threshold_p.value\": ",fitObj$threshold_p.value,",\n",sep = "")
  }else{
    cat(indent,"\t\"threshold_p.value\": null,\n",sep = "")
  }
  tableToJSON(fitObj$catalogues,tablename = "catalogues",nindent = nindent + 1)
  cat(",\n")
  tableToJSON(fitObj$signatures,tablename = "signatures",nindent = nindent + 1)
  cat(",\n")
  tableToJSON(t(fitObj$exposures),tablename = "exposures",nindent = nindent + 1)
  cat(",\n")
  vectorToJSON(fitObj$unassigned_muts,"unassigned_muts",nindent = nindent + 1)
  cat(",\n")
  vectorToJSON(fitObj$unassigned_muts_perc,"unassigned_muts_perc",nindent = nindent + 1)
  cat(",\n")
  if(all(!is.na(fitObj$bootstrap_exposures_samples))){
    cat(indent,"\t\"bootstrap_exposures_samples\": {\n",sep = "")
    for (i in 1:length(fitObj$bootstrap_exposures_samples)) {
      tableToJSON(fitObj$bootstrap_exposures_samples[[i]],tablename = rownames(fitObj$exposures)[i],nindent = nindent + 2)
      if(i<length(fitObj$bootstrap_exposures_samples)) cat(",")
      cat("\n")
    }
    cat(indent,"\t}",sep = "")
  }else{
    cat(indent,"\t\"bootstrap_exposures_samples\": null",sep = "")
  }
  cat(",\n")
  if(all(!is.na(fitObj$bootstrap_exposures_pvalues))){
    tableToJSON(fitObj$bootstrap_exposures_pvalues,tablename = "bootstrap_exposures_pvalues",nindent = nindent + 1)
  }else{
    cat(indent,"\t\"bootstrap_exposures_pvalues\": null",sep = "")
  }
  cat("\n")
  cat(indent,"}",sep = "")
  if(!disableFileWrite) {
    cat("\n")
    sink()
  }
}



#' Export the results from the FitMS function to JSON
#'
#' Exports of the results obtained with the FitMS function to JSON.
#'
#' @param fitObj object obtained from the FitMS function
#' @param filename file name where to save the JSON string
#' @param nindent number of tabs of indentation to be used, default=0. This is useful in case one wants to embed the output in a larger JSON document
#' @param disableFileWrite if TRUE the JSON lines are simply printed to screen. This option is used when fitToJSON is called from another function
#' @param blockname if specified, the first line will show the given name as the block name, which is useful if the output is embedded in a larger JSON file
#' @export
#' @examples
#' res <- FitMS(catalogues,"Breast")
#' fitMStoJSON(res,"export.json")
fitMStoJSON <- function(fitObj,
                      filename = "export.json",
                      nindent = 0,
                      disableFileWrite = FALSE,
                      blockname = NULL){
  #plot consensus exposures file with metadata
  indent <- rep("\t",nindent)
  if(!disableFileWrite) sink(file = filename)
  if(is.null(blockname)){
    cat(indent,"{\n",sep = "")
  }else{
    cat(indent,"\"",blockname,"\": {\n",sep = "")
  }
  cat(indent,"\t\"fitAlgorithm\": ",fitObj$fitAlgorithm,",\n",sep = "")
  cat(indent,"\t\"nsamples\": ",ncol(fitObj$catalogues),",\n",sep = "")
  cat(indent,"\t\"ncommonSignatures\": ",ncol(fitObj$commonSignatures),",\n",sep = "")
  cat(indent,"\t\"nrareSignatures\": ",ncol(fitObj$rareSignatures),",\n",sep = "")
  cat(indent,"\t\"nchannels\": ",nrow(fitObj$catalogues),",\n",sep = "")
  cat(indent,"\t\"method\": \"",fitObj$method,"\",\n",sep = "")
  cat(indent,"\t\"multiStepMode\": \"",fitObj$multiStepMode,"\",\n",sep = "")
  cat(indent,"\t\"maxRareSigsPerSample\": ",fitObj$maxRareSigsPerSample,",\n",sep = "")
  if(!is.null(fitObj$organ)){
    cat(indent,"\t\"organ\": \"",fitObj$organ,"\",\n",sep = "")
  }else{
    cat(indent,"\t\"organ\": null,\n",sep = "")
  }
  cat(indent,"\t\"rareSignatureTier\": \"",fitObj$rareSignatureTier,"\",\n",sep = "")
  cat(indent,"\t\"minErrorReductionPerc\": ",fitObj$minErrorReductionPerc,",\n",sep = "")
  cat(indent,"\t\"minCosSimIncrease\": ",fitObj$minCosSimIncrease,",\n",sep = "")
  cat(indent,"\t\"residualNegativeProp\": ",fitObj$residualNegativeProp,",\n",sep = "")
  cat(indent,"\t\"minCosSimRareSig\": ",fitObj$minCosSimRareSig,",\n",sep = "")

  cat(indent,"\t\"exposureFilterType\": \"",fitObj$exposureFilterType,"\",\n",sep = "")
  if(!is.na(fitObj$threshold_percent)){
    cat(indent,"\t\"threshold_percent\": ",fitObj$threshold_percent,",\n",sep = "")
  }else{
    cat(indent,"\t\"threshold_percent\": null,\n",sep = "")
  }
  if(!is.na(fitObj$giniThresholdScaling)){
    cat(indent,"\t\"giniThresholdScaling\": ",fitObj$giniThresholdScaling,",\n",sep = "")
  }else{
    cat(indent,"\t\"giniThresholdScaling\": null,\n",sep = "")
  }
  cat(indent,"\t\"useBootstrap\": ",ifelse(fitObj$useBootstrap,"true","false"),",\n",sep = "")
  if(!is.na(fitObj$nboot)){
    cat(indent,"\t\"nboot\": ",fitObj$nboot,",\n",sep = "")
  }else{
    cat(indent,"\t\"nboot\": null,\n",sep = "")
  }
  if(!is.na(fitObj$threshold_p.value)){
    cat(indent,"\t\"threshold_p.value\": ",fitObj$threshold_p.value,",\n",sep = "")
  }else{
    cat(indent,"\t\"threshold_p.value\": null,\n",sep = "")
  }
  if(length(fitObj$whichSamplesMayHaveRareSigs)>0){
    vectorToJSON(fitObj$whichSamplesMayHaveRareSigs,"whichSamplesMayHaveRareSigs",nindent = nindent + 1)
    cat(",\n")
  }else{
    cat(indent,"\t\"whichSamplesMayHaveRareSigs\": null,\n",sep = "")
  }
  if(length(fitObj$rareSigChoice)>0){
    cat(indent,"\t\"rareSigChoice\": {\n",sep = "")
    for (i in 1:length(fitObj$rareSigChoice)){
      cat(indent,"\t\t\"",names(fitObj$rareSigChoice)[i],"\": \"",fitObj$rareSigChoice[[names(fitObj$rareSigChoice)[i]]],"\"",sep = "")
      if(i < length(fitObj$rareSigChoice)) cat(",")
      cat("\n")
    }
    cat(indent,"\t},\n",sep = "")
  }else{
    cat(indent,"\t\"rareSigChoice\": null,\n",sep = "")
  }
  if(length(fitObj$candidateRareSigsCosSim)>0){
    cat(indent,"\t\"candidateRareSigsCosSim\": {\n",sep = "")
    for (i in 1:length(fitObj$candidateRareSigsCosSim)){
      sampleName <- names(fitObj$candidateRareSigsCosSim)[i]
      candidates <- fitObj$candidateRareSigsCosSim[[sampleName]]
      cat(indent,"\t\t\"",sampleName,"\": {\n",sep = "")
      for (j in 1:length(candidates)) {
        cat(indent,"\t\t\t\"",names(candidates)[j],"\": ",candidates[j],sep = "")
        if(j < length(fitObj$candidates)) cat(",")
        cat("\n")
      }
      cat(indent,"\t\t}",sep = "")
      if(i < length(fitObj$rareSigChoice)) cat(",")
      cat("\n")
    }
    cat(indent,"\t},\n",sep = "")
  }else{
    cat(indent,"\t\"candidateRareSigsCosSim\": null,\n",sep = "")
  }
  tableToJSON(fitObj$catalogues,tablename = "catalogues",nindent = nindent + 1)
  cat(",\n")
  tableToJSON(fitObj$commonSignatures,tablename = "commonSignatures",nindent = nindent + 1)
  cat(",\n")
  tableToJSON(fitObj$rareSignatures,tablename = "rareSignatures",nindent = nindent + 1)
  cat(",\n")
  tableToJSON(t(fitObj$exposures),tablename = "exposures",nindent = nindent + 1)
  cat(",\n")
  if(all(!is.na(fitObj$bootstrap_exposures_samples))){
    cat(indent,"\t\"bootstrap_exposures_samples\": {\n",sep = "")
    for (i in 1:length(fitObj$bootstrap_exposures_samples)) {
      tableToJSON(fitObj$bootstrap_exposures_samples[[i]],tablename = names(fitObj$bootstrap_exposures_samples)[i],nindent = nindent + 2)
      if(i<length(fitObj$bootstrap_exposures_samples)) cat(",")
      cat("\n")
    }
    cat(indent,"\t}",sep = "")
  }else{
    cat(indent,"\t\"bootstrap_exposures_samples\": \"NULL\"",sep = "")
  }
  cat(",\n")

  cat(indent,"\t\"all_fits_samples\": {\n",sep = "")
  for (i in 1:length(fitObj$samples)) {
    sampleName <- names(fitObj$samples)[i]
    contentNames <- names(fitObj$samples[[sampleName]])
    cat(indent,"\t\t\"",sampleName,"\": {\n",sep = "")
    fitToJSON(fitObj$samples[[sampleName]]$fitCommonOnly,disableFileWrite = T,blockname = "fitCommonOnly",nindent = nindent + 3)
    if("fitWithRare" %in% contentNames){
      cat(",\n")
      rareNames <- names(fitObj$samples[[sampleName]][["fitWithRare"]])
      cat(indent,"\t\t\t\"fitWithRare\": {\n",sep = "")
      for(j in 1:length(rareNames)){
        fitToJSON(fitObj$samples[[sampleName]]$fitWithRare[[rareNames[j]]],disableFileWrite = T,blockname = rareNames[j],nindent = nindent + 4)
        if(j<length(rareNames)) cat(",")
        cat("\n")
      }
      cat(indent,"\t\t\t}",sep = "")
    }
    cat("\n")
    cat(indent,"\t\t}",sep = "")
    if(i<length(fitObj$samples)) cat(",")
    cat("\n")
  }
  cat(indent,"\t}",sep = "")

  cat("\n")
  cat(indent,"}\n",sep = "")
  if(!disableFileWrite) sink()
}

#' Save fit object to file
#'
#' This function saves a Fit or FitMS object to an R data file using a standard name
#'
#' @param fitObj object obtained from the Fit of FitMS function
#' @param filename file name where to save the fit object
#' @export
#' @examples
#' fitObj <- FitMS(catalogues,organ="Breast")
#' saveFitToFile(fitObj,"fit.rData")
saveFitToFile <- function(fitObj,filename,verbose=T){
  if(verbose) message("[info saveFitToFile] saving ",fitObj$fitAlgorithm," object to file ",filename)
  save(file = filename,fitObj)
}

#' Load fit object from file
#'
#' This function loads a Fit or FitMS object from an R data file.
#' In order to work, the file must have been generated using the saveFitToFile function,
#' which will save the fit object with a standard name.
#'
#' @param filename file name from which to load the fit object
#' @export
#' @examples
#' fitObj <- loadFitFromFile("fit.rData")
loadFitFromFile <- function(filename,verbose=T){
  load(file = filename)
  if(verbose) message("[info loadFitToFile] loaded ",fitObj$fitAlgorithm," object from file ",filename)
  return(fitObj)
}
