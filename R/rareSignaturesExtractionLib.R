#' Rare Signatures Extraction
#'
#' Extract rare signatures from sample catalogues. This function is the last step of
#' the extraction of rare signatures from sample catalogues. Before running this
#' function one should use the unexplainedSamples function to identify which
#' samples are likely to contain a rare signature. The unexplainedSamples function
#' produces a residual for each sample, indicating the part of the catalogue that is
#' not explained by the common signatures used. Then, one should cluster the residuals
#' of the samples that the unexplainedSamples function considered significant, using the
#' cataloguesClustering function. Finally, one should be ready to run this function.
#' See Examples below. This function will use the NNLM package to extract rare functions
#' while at the same time fitting the common signatures. An extraction is performed
#' for each of the clusters of residuals obtained using cataloguesClustering.
#' The extraction will use a variant of the Lee and Seung multiplicative algorithm,
#' optimising the KLD. The starting point for the extraction will be the commonExposures
#' provided. Specific signature exposures can be ignored with the parameter commonSigsToIgnore.
#'
#' @param outfileRoot if specified, generate a plot, otherwise no plot is generated
#' @param catalogues original catalogues, channels as rows and samples as columns
#' @param commonSignatures common mutational signatures used for fitting, channels as rows, signatures as columns
#' @param commonExposures exposures obtained from fitting catalogues using commonSignatures, signatures names as row names and sample names as column names. Typically this should be obtained from the unexplainedSamples function.
#' @param residuals residuals from catalogues fitted using common signatures. Typically this should be obtained from the unexplainedSamples function.
#' @param unexpl_samples names of samples that were found to be unexplained by the unexplainedSamples function
#' @param clusters vector of integers indicating the cluster that each unexplained sample belongs to. Typically this is obtained clustering residuals of unexplained samples using the catalogueClustering function.
#' @param useclusters list object indicating which cluster or group of clusters should be used to extract each rare signature. This is a list object and the length of the list indicates how many rare signatures will be extracted. Each entry of the list is a vector indicating one or more cluster numbers that will be combined into one for a rare signature extraction. For example, list(c(1,4),c(3),c(2)) indicates that three rare signatures will be extracted, one using samples from clusters 1 and 4, one using samples in cluster 3, and one using samples in cluster 2.
#' @param commonSigsToIgnore this parameter is used to reduce the importance of specific common signatures, by setting their exposures to 0 in the commonExposures table, thus changing the starting point of the optimisation. This is very useful when the rare signature to be extracted bears similarities with common signatures, which in turn interfere with the extraction. Leave NULL if not used, or use NA for extraction where this should not be used. For example, list(NA,c("commonSig1","commonSig3"),NA) assumes that three rare signatures will be extracted (determined by the useclusters parameter) and that the second extraction should ignore commonSig1 and commonSig3, which should also be row names of commonExposures.
#' @param checkForMissed if set to TRUE, there will be an additional check where each extracted rare signature will be compared to all the provided residuals to identify additional samples that may have the rare signature (min cosine similarity signature vs residual 0.95)
#' @param maxiter maximum number of iterations for the Lee and Seung multiplicative algorithm. Can be a single value or a vector if different extractions should have different maxiter values. It should be noted that performing an extraction on a single sample has potentially infinite solutions. For this reason it is useful in this case to run the Lee and Seung multiplicative algorithm starting from the commonExposures solution and run a limited number of iterations to find a neighbouring solution
#' @param nsigs How many signatures should be extracted. Typically this should be left to 1, indicating that for each group of unexplained samples (indicated by each entry in useclusters), only one rare signature should be extracted. If multiple rare signatures are thought to be present in one of the groups, then nsigs can be a vector (same length as useclusters) and a different values can be used, for example c(1,2,1) will extract 2 rare signatures from the samples in the second group.
#' @keywords extract rare signatures
#' @return list of rare signatures and list of sample names for each signature
#' @export
#' @examples
#' resUnexpl <- unexplainedSamples(catalogues=catalogues,
#'                                 sigs=signatures)
#' significant_residuals <- resUnexpl$all_residuals[,resUnexpl$which_significant]
#' clusterResiduals <- cataloguesClustering(significant_residuals,
#'                                          nclusters = 1:5)
#' resObj <- rareSignatureExtraction(catalogues=catalogues,
#'                                   commonSignatures=signatures,
#'                                   commonExposures=resUnexpl$exposures,
#'                                   residuals=resUnexpl$all_residuals,
#'                                   unexpl_samples=resUnexpl$unexplSamples,
#'                                   clusters=clusterResiduals$clusters_table[,"3"],
#'                                   useclusters=list=(c(1),c(3)),
#'                                   maxiter=c(100,1000))
rareSignatureExtraction <- function(outfileRoot,
                                    catalogues,
                                    commonSignatures,
                                    commonExposures,
                                    residuals,
                                    unexpl_samples,
                                    clusters,
                                    useclusters,
                                    commonSigsToIgnore=NULL,
                                    checkForMissed=FALSE,
                                    maxiter=1000,
                                    nsigs=1){
  listofsignatures <- list()
  listofsamples <- list()

  for (rarecounter in 1:length(useclusters)) {
    message("[info rareSignatureExtraction] Extracting rare signatures: ", rarecounter, " of ",length(useclusters))

    currentSampleNames <- unexpl_samples[clusters %in% useclusters[[rarecounter]]]
    currentSamples <- catalogues[,currentSampleNames,drop=F]
    currentResiduals <- residuals[,currentSampleNames,drop=F]
    expoCurrentSamples <- commonExposures[,currentSampleNames,drop=F]

    # set common sigs to ignore if any
    if(typeof(commonSigsToIgnore)=="list"){
      if(!is.null(commonSigsToIgnore[[rarecounter]]) & !is.na(commonSigsToIgnore[[rarecounter]])){
        expoCurrentSamples[commonSigsToIgnore[[rarecounter]],] <- 0
      }
    }
    # check max iter
    if(length(maxiter)==1){
      usemaxiter <- maxiter
    }else{
      usemaxiter <- maxiter[rarecounter]
    }
    # check nsigs
    if(length(nsigs)==1){
      usensigs <- nsigs
    }else{
      usensigs <- nsigs[rarecounter]
    }
    suppressWarnings(
      res.nnmf <- NNLM::nnmf(as.matrix(currentSamples),
                             init = list(W0 = as.matrix(commonSignatures),H1 = as.matrix(expoCurrentSamples)),
                             loss = "mkl",
                             method = "lee",
                             k = usensigs,
                             max.iter = usemaxiter,
                             check.k = FALSE)
    )
    rareSignatures <- apply(res.nnmf$W[,1:usensigs,drop=F],2,function(x)x/sum(x))
    signames <- paste0("U",rarecounter)
    if(usensigs>1) signames <- paste0(signames,letters[1:usensigs])
    colnames(rareSignatures) <- signames

    #listofsignatures[[paste0("rare",rarecounter)]] <- rareSignatures
    #listofsamples[[paste0("rare",rarecounter)]] <- colnames(currentSamples)

    samplesMissed <- c()
    if(checkForMissed){
      residualPositiveUnseen <- residuals[,setdiff(colnames(residuals),unexpl_samples),drop=F]
      residualPositiveUnseen[residualPositiveUnseen<0] <- 0
      resCorr <- computeCorrelationOfTwoSetsOfSigs(residualPositiveUnseen,rareSignatures)
      whichSimilar <- rownames(resCorr)[which(resCorr>=0.95)]
      if(length(whichSimilar)>0){
        #plotDNVSignatures(R[,whichSimilar,drop=F])
        #plotDNVSignatures(catalogue_full[,whichSimilar,drop=F])
        samplesMissed <- whichSimilar
        message("[info rareSignatureExtraction] ",signames," recovered missed samples: ",paste(whichSimilar,collapse = ", "))
      }
    }
    currentSampleNames <- union(currentSampleNames,samplesMissed)
    listofsignatures[[paste0("rare",rarecounter)]] <- rareSignatures
    listofsamples[[paste0("rare",rarecounter)]] <- currentSampleNames

    # update
    currentSamples <- catalogues[,currentSampleNames,drop=F]
    currentResiduals <- residuals[,currentSampleNames,drop=F]
  }

  # combine to save and plot
  fullorgansigs <- cbind(commonSignatures,do.call(cbind,listofsignatures))

  # plotting
  if(!is.null(outfileRoot)){



    # create the directory if missing
    trimdir <- function(x){
      while(substr(x,nchar(x),nchar(x))!="/" & nchar(x)>0){
        x <- substr(x,1,nchar(x)-1)
      }
      return(x)
    }

    # create dir
    xdir <- trimdir(outfileRoot)
    if(nchar(xdir)>0) dir.create(xdir,recursive = T,showWarnings = F)

    message("[info rareSignatureExtraction] Saving to file in ",ifelse(nchar(xdir)>0,xdir,"current directory"),".")

    writeTable(fullorgansigs[,(ncol(commonSignatures) + 1):(ncol(fullorgansigs)),drop=F],paste0(outfileRoot,"_","signatures_rare.tsv"))
    plotSignatures(fullorgansigs[,(ncol(commonSignatures) + 1):(ncol(fullorgansigs)),drop=F],paste0(outfileRoot,"_","signatures_rare.pdf"),plot_sum = F,ncolumns = 3)
    # saving
    writeTable(fullorgansigs,paste0(outfileRoot,"_","signatures.tsv"))
    plotSignatures(fullorgansigs,paste0(outfileRoot,"_","signatures.pdf"),plot_sum = F,ncolumns = 3)

    for (rarecounter in 1:length(useclusters)) {
      currentSampleNames <- listofsamples[[paste0("rare",rarecounter)]]
      currentSamples <- catalogues[,currentSampleNames,drop=F]
      currentResiduals <- residuals[,currentSampleNames,drop=F]
      rareSignatures <- listofsignatures[[paste0("rare",rarecounter)]]
      plotSignatures(currentSamples,output_file = paste0(outfileRoot,"_","SamplesCluster",rarecounter,"_Catalogues.pdf"),ncolumns = 3)
      plotSignatures(currentResiduals,output_file = paste0(outfileRoot,"_","SamplesCluster",rarecounter,"_Residuals.pdf"),ncolumns = 3)
      plotSignatures(rareSignatures,output_file = paste0(outfileRoot,"_","SamplesCluster",rarecounter,"_rareSignatures.pdf"),plot_sum = F,ncolumns = 3)
    }
  }

  res <- list()
  res$listofsignatures <- listofsignatures
  res$listofsamples <- listofsamples
  res$commonAndRareSignatures <- fullorgansigs
  return(res)
}



#' Finalise Common and Rare Signature Exposures
#'
#' This function is part of the pipeline for the extraction of rare signatures. After the
#' rare signatures have been extracted and the corresponding samples identified
#' using the rareSignatureExtraction function, this function fits the common and
#' rare signatures into the sample catalogues and reports also a comparison
#' between the normalised error withvs without rare signatures.
#'
#'
#' @param outfileRoot output directory and base file name, for example outdir/BaseName. Multiple output files will be generated with filenames the chosen BaseName as prefix.
#' @param catalogues original catalogues, channels as rows and samples as columns
#' @param commonSigs common mutational signatures used for fitting, channels as rows, signatures as columns
#' @param listofsignatures listofsignatures object obtained from the rareSignatureExtraction function. It contains the rare signatures extracted
#' @param listofsamples listofsamples object obtained from the rareSignatureExtraction function. It contains the samples that have each of the rare signatures extracted
#' @param nboot number of bootstraps for the signature fit
#' @param threshold_percent minimum percentage of mutations in a sample that should be assigned to a signature
#' @param threshold_p.value p-value to decide whether a signature exposure is statistically higher then the threshold_percent
#' @param min_sample_muts minimum number of mutations in a sample to perform signature analysis. If a sample has less than min_sample_muts mutations, then all signature exposures are set to zero and all mutations are unsassigned
#' @param nparallel number of parallel processes to use for the signature fit
#' @keywords fit rare common signatures
#' @return exposures and errors for the signature fit using common+rare signatures and common only signatures
#' @export
#' @examples
#' resUnexpl <- unexplainedSamples(catalogues=catalogues,
#'                                 sigs=signatures)
#' significant_residuals <- resUnexpl$all_residuals[,resUnexpl$which_significant]
#' clusterResiduals <- cataloguesClustering(significant_residuals,
#'                                          nclusters = 1:5)
#' resObj <- rareSignatureExtraction(catalogues=catalogues,
#'                                   commonSignatures=signatures,
#'                                   commonExposures=resUnexpl$exposures,
#'                                   residuals=resUnexpl$all_residuals,
#'                                   unexpl_samples=resUnexpl$unexplSamples,
#'                                   clusters=clusterResiduals$clusters_table[,"3"],
#'                                   useclusters=list=(c(1),c(3)),
#'                                   maxiter=c(100,1000))
#' finaliseCommonRareSignatureExposures(outfileRoot = "outdir/BaseName",
#'                                      catalogues = catalogues,
#'                                      commonSigs = signatures,
#'                                      listofsignatures = resObj$listofsignatures,
#'                                      listofsamples = resObj$listofsamples)
finaliseCommonRareSignatureExposures <- function(outfileRoot,
                                                 catalogues,
                                                 commonSigs,
                                                 listofsignatures=NULL,
                                                 listofsamples=NULL,
                                                 nboot = 200,
                                                 threshold_percent = 5,
                                                 threshold_p.value = 0.05,
                                                 min_sample_muts = 100,
                                                 nparallel = 1){

  if(is.null(outfileRoot)){
    message("[error finaliseCommonRareSignatureExposures] please provide the output directy and root file name using the outfileRoot parameter. Nothing done.")
    return(NULL)
  }

  # combine signatures
  fullorgansigs <- commonSigs
  if(!is.null(listofsignatures)) fullorgansigs <- cbind(fullorgansigs,do.call(cbind,listofsignatures))

  message("[info finaliseCommonRareSignatureExposures] setting up the fit mask.")

  # prepare the signature fit mask
  sigfitmask <- matrix(TRUE,ncol = ncol(catalogues),nrow = ncol(commonSigs),dimnames = list(colnames(commonSigs),colnames(catalogues)))
  sigfitmaskCommonOnly <- matrix(TRUE,ncol = ncol(catalogues),nrow = ncol(commonSigs),dimnames = list(colnames(commonSigs),colnames(catalogues)))

  if(!is.null(listofsignatures)){
    newsigcount <- ncol(commonSigs) + 1

    for (i in names(listofsignatures)){
      # i <- "rare1"
      newsigs <- listofsignatures[[i]]
      newsigssamples <- listofsamples[[i]]

      newsigcount <- newsigcount + ncol(newsigs)

      newrows <- matrix(FALSE,nrow = ncol(newsigs),ncol = ncol(sigfitmask),dimnames = list(colnames(newsigs),colnames(sigfitmask)))
      # sigfitmaskCommonOnly <- rbind(sigfitmaskCommonOnly,newrows)
      newrows[,newsigssamples] <- TRUE
      sigfitmask <- rbind(sigfitmask,newrows)
    }
  }

  message("[info finaliseCommonRareSignatureExposures] Fitting using common and rare signatures...")

  resfitWithRare <- fitSignaturesAndSave(outfileRoot = paste0(outfileRoot,"_fitWithRare"),
                                         fullorgansigs = fullorgansigs,
                                         sigfitmask = sigfitmask,
                                         catalogues = catalogues,
                                         nboot = nboot,
                                         threshold_percent = threshold_percent,
                                         threshold_p.value = threshold_p.value,
                                         min_sample_muts = min_sample_muts,
                                         nparallel = nparallel)

  message("[info finaliseCommonRareSignatureExposures] Fitting using common signatures only...")

  resfitCommonOnly <- fitSignaturesAndSave(outfileRoot = paste0(outfileRoot,"_fitCommonOnly"),
                                           fullorgansigs = commonSigs,
                                           sigfitmask = sigfitmaskCommonOnly,
                                           catalogues = catalogues,
                                           nboot = nboot,
                                           threshold_percent = threshold_percent,
                                           threshold_p.value = threshold_p.value,
                                           min_sample_muts = min_sample_muts,
                                           nparallel = nparallel)


  message("[info finaliseCommonRareSignatureExposures] Comparing common only with common+rare...")

  # compare results
  resCompare <- compareError(outfileRoot = outfileRoot,
                             catalogues = catalogues,
                             sigs1 = commonSigs,
                             exposures1 = resfitCommonOnly$exposures[-nrow(resfitCommonOnly$exposures),colnames(catalogues)],
                             sigs2 = fullorgansigs,
                             exposures2 = resfitWithRare$exposures[-nrow(resfitWithRare$exposures),colnames(catalogues)],
                             samplesWithRareSig = unlist(listofsamples))

  # plotting

  # create the directory if missing
  trimdir <- function(x){
    while(substr(x,nchar(x),nchar(x))!="/" & nchar(x)>0){
      x <- substr(x,1,nchar(x)-1)
    }
    return(x)
  }

  # this is where the plots/tables output go
  writeTable(sigfitmask,paste0(outfileRoot,"_","sigfitmask.tsv"))
  plotExposures(resfitCommonOnly$exposures,output_file = paste0(outfileRoot,"_fitCommonOnlyExp.pdf"))
  plotExposures(resfitWithRare$exposures,output_file = paste0(outfileRoot,"_fitWithRareExp.pdf"))

  message("[info finaliseCommonRareSignatureExposures] done.")

  # any return obj
  res <- list()
  res$sigfitmask <- sigfitmask
  res$sigfitmaskCommonOnly <- sigfitmaskCommonOnly
  res$fitWithRare <- resfitWithRare
  res$fitCommonOnly <- resfitCommonOnly
  res$compareErrorCommonOnlyVsCommonAndRare <- resCompare
  return(res)
}




fitSignaturesAndSave <- function(outfileRoot,
                                 fullorgansigs,
                                 sigfitmask,
                                 catalogues,
                                 nboot = 200,
                                 threshold_percent = 5,
                                 threshold_p.value = 0.05,
                                 min_sample_muts = 100,
                                 nparallel = 1){
  res.fit.boot <- matrix(0,ncol = ncol(sigfitmask),nrow = nrow(sigfitmask),dimnames = dimnames(sigfitmask))
  # add unassigned mutations row
  res.fit.boot <- rbind(res.fit.boot,matrix(0,nrow = 1,ncol = ncol(res.fit.boot),dimnames = list("unassigned",colnames(res.fit.boot))))

  res.boot.data <- list()

  for (i in 1:ncol(sigfitmask)){
    # i <- 1
    message("sample ",i, " of ",ncol(sigfitmask))
    sigstouse <- fullorgansigs[,sigfitmask[,i],drop=F]
    sampleMuts <- sum(catalogues[,i,drop=F])
    if(sampleMuts>=min_sample_muts){
      res.boot.data[[i]] <- SignatureFit_withBootstrap(cat = catalogues[,i,drop=F],
                                                                            signature_data_matrix = sigstouse,
                                                                            nboot = nboot,
                                                                            threshold_percent = threshold_percent,
                                                                            threshold_p.value = threshold_p.value,
                                                                            verbose = F,
                                                                            doRound = F,
                                                                            showDeprecated = F,
                                                                            nparallel = nparallel)
      res.fit.boot[rownames(res.boot.data[[i]]$E_median_filtered),i] <- res.boot.data[[i]]$E_median_filtered
      # add unassigned
      res.fit.boot["unassigned",i] <- res.boot.data[[i]]$unassigned_muts
    }else{
      res.fit.boot[,i] <- 0
      res.fit.boot["unassigned",i] <- sampleMuts
    }
  }

  writeTable(res.fit.boot,file = paste0(outfileRoot,"_","exposures_final.tsv"))
  save(file = paste0(outfileRoot,"_","exposures_bootstraps.rData"),res.boot.data)

  res <- list()
  res$exposures <- res.fit.boot
  res$bootstraps <- res.boot.data
  return(res)
}


compareError <- function(outfileRoot=NULL,
                         catalogues,
                         sigs1,
                         exposures1,
                         sigs2,
                         exposures2,
                         samplesWithRareSig){
  error1 <- catalogues - as.matrix(sigs1) %*% as.matrix(exposures1)
  error2 <- catalogues - as.matrix(sigs2) %*% as.matrix(exposures2)
  errorSAD1 <- c()
  errorSAD2 <- c()
  norm_errorSAD1 <- c()
  norm_errorSAD2 <- c()
  for (i in 1:ncol(catalogues)){
    errorSAD1 <- c(errorSAD1,sum(abs(error1[,i])))
    errorSAD2 <- c(errorSAD2,sum(abs(error2[,i])))
    sumcati <- sum(catalogues[,i])
    norm_errorSAD1 <- c(norm_errorSAD1,errorSAD1[i]/sumcati)
    norm_errorSAD2 <- c(norm_errorSAD2,errorSAD2[i]/sumcati)
  }
  plottitle <- "Normalised Sum of Absolute Deviations (SAD)\nbetween original and reconstructed catalogue"
  xlab <- "SAD/n. mutations (common+rare sigs.)"
  ylab <- "SAD/n. mutations (common sigs. only)"
  cairo_pdf(filename = paste0(outfileRoot,"_compareErrorAcrossSamples.pdf"),width = 6,height = 5,pointsize = 14)
  par(mar=c(4,6,3.5,2),mgp=c(2.5,1,0))
  plot(norm_errorSAD2,norm_errorSAD1,
       xlim = c(0,max(norm_errorSAD2)),
       ylim = c(0,max(norm_errorSAD1)),
       col=rgb(0.5,0.5,0.5,0.5),
       pch = 16,asp=1,
       xlab = xlab,las=1,
       ylab = ylab,
       main = plottitle,cex.main=0.9)
  points(norm_errorSAD2[colnames(catalogues) %in% samplesWithRareSig],
         norm_errorSAD1[colnames(catalogues) %in% samplesWithRareSig],
         col=rgb(0.6,0,0,1),pch = 16)
  legend(x="bottomright",col = rgb(0.6,0,0,1),pch = 16,legend = "samples with rare sig",bty = 'n',cex = 0.9)
  abline(a = 0, b = 1,lty = 2)
  dev.off()

  # return object
  res <- list()
  res$norm_errorSAD1 <- norm_errorSAD1
  res$norm_errorSAD2 <- norm_errorSAD2
  return(res)
}
