

evaluatePerformanceExposuresCore <- function(true_exposures,
                                         estimated_exposures){
  true_exposures <- as.data.frame(true_exposures,stringsAsFactors = F)
  estimated_exposures <- as.data.frame(estimated_exposures,stringsAsFactors = F)
  # first fill in some blanks
  allSigNames <- union(colnames(true_exposures),colnames(estimated_exposures))
  # fix missing in original
  true_exposures[,setdiff(allSigNames,colnames(true_exposures))] <- 0
  # fix missing in the estimated
  estimated_exposures[,setdiff(allSigNames,colnames(estimated_exposures))] <- 0

  # now reorder the columns and rows
  estimated_exposures <- estimated_exposures[rownames(true_exposures),colnames(true_exposures)]

  TP <- sum(estimated_exposures[true_exposures>0]>0)
  FP <- sum(estimated_exposures[true_exposures==0]>0)
  TN <- sum(estimated_exposures[true_exposures==0]==0)
  FN <- sum(estimated_exposures[true_exposures>0]==0)

  # now compute the TP,TN,FP,FN assignments
  resPerf <- list()
  # npos <- sum(estimated_exposures>0)
  resPerf$TP <- sum(estimated_exposures[true_exposures>0]>0)
  resPerf$FP <- sum(estimated_exposures[true_exposures==0]>0)
  # nneg <- sum(estimated_exposures==0)
  resPerf$TN <- sum(estimated_exposures[true_exposures==0]==0)
  resPerf$FN <- sum(estimated_exposures[true_exposures>0]==0)
  resPerf$sensitivity <- TP/(TP+FN)
  resPerf$specificity <- TN/(TN+FP)
  resPerf$F1score <- 2*TP/(2*TP+FP+FN)
  resPerf$PPV <- TP/(TP+FP)
  resPerf$NPV <- TN/(TN+FN)
  resPerf$rmse <- sqrt(sum((true_exposures-estimated_exposures)^2)/(ncol(true_exposures)*nrow(true_exposures)))
  return(resPerf)
}

#' evaluate performance exposures
#'
#' Compare two exposure matrices, one with the true exposures and one with estimated exposures, and
#' compute various metrics to determine how well the estimated exposures match the true exposures.
#' If the names of common and/or rare mutational signatures are provided, then the
#' metrics are calculated separately for common and rare signatures.
#' Results are plotted to a file if outfile is specified.
#'
#' @param true_exposures matrix of true exposures, for example obtained using simulations, with mutational signatures names as colnames and sample names as rownames
#' @param estimated_exposures matrix of estimated exposures, for example obtained using a signature extraction or fit algorithm on a simulated dataset, with mutational signatures names as colnames and sample names as rownames
#' @param commonNames names of common mutational signatures. If NULL, all signatures in the colnames of true_exposures and estimated_exposures will be used
#' @param rareNames names of rare mutational signatures. If NULL, rare signatures metrics will not be computed. This is ok if the simulated data didn't have any rare signatures.
#' @param outfile file name for plotting, please use .pdf file name extension. Can be omitted.
#' @return performance metrics
#' @export
evaluatePerformanceExposures <- function(true_exposures,
                                         estimated_exposures,
                                         commonNames = NULL,
                                         rareNames = NULL,
                                         outfile = NULL){
  if(!is.null(commonNames)){
    resPerfCommon <- evaluatePerformanceExposuresCore(true_exposures[,intersect(commonNames,colnames(true_exposures)),drop=F],
                                                      estimated_exposures[,intersect(commonNames,colnames(estimated_exposures)),drop=F])
  }else{
    # the core function will take care of merging the signature names
    resPerfCommon <- evaluatePerformanceExposuresCore(true_exposures,
                                                      estimated_exposures)
  }

  if(!is.null(rareNames)){
    resPerfRare <- evaluatePerformanceExposuresCore(true_exposures[,intersect(rareNames,colnames(true_exposures)),drop=F],
                                                    estimated_exposures[,intersect(rareNames,colnames(estimated_exposures)),drop=F])
  }else{
    resPerfRare <- NULL
  }

  perfList <- list()
  perfList[["performance"]] <- list(resPerfCommon=resPerfCommon,
                                    resPerfRare=resPerfRare)
  # convert to tables
  perfTables <- combinePerformanceExposures(perfList)

  if(!is.null(outfile)){
    plotPerformanceExposures(perfTables,outfile)
  }

  return(perfTables)
}

#' evaluate performance exposures list
#'
#' Compare an exposure matrix with the true exposures with a list of estimated exposure matrices, and
#' compute various metrics to determine how well the estimated exposures match the true exposures.
#' If the names of common and/or rare mutational signatures are provided, then the
#' metrics are calculated separately for common and rare signatures.
#' Results are plotted to a file if outfile is specified.
#'
#' @param true_exposures matrix of true exposures, for example obtained using simulations, with mutational signatures names as colnames and sample names as rownames
#' @param estimated_exposures_list list of matrices of estimated exposures, for example obtained using a signature extraction or fit algorithm on a simulated dataset, with mutational signatures names as colnames and sample names as rownames. Use a named list, so that the names can be used in he data structure returned and plots.
#' @param commonNames names of common mutational signatures. If NULL, all signatures in the colnames of true_exposures and estimated_exposures will be used
#' @param rareNames names of rare mutational signatures. If NULL, rare signatures metrics will not be computed. This is ok if the simulated data didn't have any rare signatures.
#' @param outfile file name for plotting, please use .pdf file name extension. Can be omitted.
#' @return performance metrics
#' @export
evaluatePerformanceExposuresList <- function(true_exposures,
                                             estimated_exposures_list,
                                             commonNames = NULL,
                                             rareNames = NULL,
                                             outfile = NULL){
  groupNames <- names(estimated_exposures_list)
  perfList <- list()
  for (g in groupNames) {
    perfList[[g]] <- evaluatePerformanceExposures(true_exposures,
                                                  estimated_exposures_list[[g]],
                                                  commonNames = commonNames,
                                                  rareNames = rareNames)
  }

  # combine tables
  perfTables <- list()
  for (g in groupNames) {
    for (n in names(perfList[[g]])){
      if(ncol(perfList[[g]][[n]])>0){
        colnames(perfList[[g]][[n]]) <- g
        if(is.null(perfTables[[n]])){
          perfTables[[n]] <- perfList[[g]][[n]]
        }else{
          perfTables[[n]] <- cbind(perfTables[[n]],perfList[[g]][[n]])
        }
      }
    }
  }

  if(!is.null(outfile)){
    plotPerformanceExposures(perfTables,outfile)
  }

  return(perfTables)
}


combinePerformanceExposures <- function(perfList){
  perfNames <- c()
  perfGroups <- c()
  for (n in names(perfList)){
    perfGroups <- union(perfGroups,names(perfList[[n]]))
    for (ni in names(perfList[[n]])) {
      perfNames <- union(perfNames,names(perfList[[n]][[ni]]))
    }
  }
  returnTables <- list()
  for(g in perfGroups) returnTables[[g]] <- data.frame(row.names = perfNames,check.names = F,stringsAsFactors = F)
  for (n in names(perfList)){
    for (ni in names(perfList[[n]])) {
      for(nii in names(perfList[[n]][[ni]])) returnTables[[ni]][nii,n] <- perfList[[n]][[ni]][[nii]]
    }
  }
  return(returnTables)
}


plotPerformanceExposures <- function(perfTables,
                                     outfile = NULL){

  dir.create(dirname(outfile),showWarnings = F,recursive = T)

  nmetrics <- nrow(perfTables$resPerfCommon)
  if(is.null(perfTables$resPerfRare)) {
    plotRareSigs <- FALSE
  }else if (!(ncol(perfTables$resPerfRare)==0)){
    plotRareSigs <- TRUE
    # also replace NA and NaN with zeros
    perfTables$resPerfRare[is.na(perfTables$resPerfRare)] <- 0
  }
  allcolours <- c("#F3C300","#875692","#F38400","#A1CAF1","#BE0032","#C2B280","#848482","#008856","#E68FAC","#0067A5","#F99379","#604E97","#F6A600","#B3446C","#DCD300")

  pointsize <- 16
  maxNamesLength <- max(sapply(colnames(perfTables$resPerfCommon),function(x) strwidth(x,units = "inch",ps = pointsize),USE.NAMES = F))

  plotHeight <- 5*nmetrics
  plotWidth <- ifelse(plotRareSigs,12,6)

  cairo_pdf(filename = outfile,height = plotHeight,width = plotWidth,pointsize = pointsize)
  par(mfrow=c(nmetrics,ifelse(plotRareSigs,2,1)),mai=c(maxNamesLength+0.2,1,1,1))
  for(i in 1:nmetrics){
    metricsName <- rownames(perfTables$resPerfCommon)[i]
    if(plotRareSigs) metricsName <- paste0(metricsName," - common")
    dataToPlot <- unlist(perfTables$resPerfCommon[i,,drop=T])
    ytop <- max(dataToPlot)
    if(ytop<=1) ytop <- 1
    barplot(dataToPlot,beside = T,ylim = c(0,ytop),
            names.arg = colnames(perfTables$resPerfCommon),
            col = allcolours[1:ncol(perfTables$resPerfCommon)],border = NA,
            main = metricsName,las=2)
    if(plotRareSigs) {
      dataToPlot <- unlist(perfTables$resPerfRare[i,,drop=T])
      ytop <- max(dataToPlot)
      if(ytop<=1) ytop <- 1
      barplot(dataToPlot,beside = T,ylim = c(0,ytop),
              names.arg = colnames(perfTables$resPerfRare),
              col = allcolours[1:ncol(perfTables$resPerfRare)],border = NA,
              main = paste0(rownames(perfTables$resPerfRare)[i]," - rare"),las=2)
    }
  }
  dev.off()
  return(perfTables)
}


#' evaluate performance signature similarity
#'
#' Compare a signatures matrix of true signatures with a matrix of estimated signatures, and
#' compute a match between the signatures and their cosine similarity. The match implies
#' an optimal assignment that maximises the cosine similarity while allowing each true_signature
#' to match at most one estimated signature and viceversa.
#' If the number of signatures differ, the least similar signatures are removed until the number
#' of signatures is the same
#' Results are plotted to a file if outfile is specified.
#'
#' @param true_signatures matrix of true signatures, for example obtained using simulations, with mutational signatures names as colnames and channels as rownames
#' @param estimated_signatures matrix of estimated signatures, for example obtained using a signature extraction on a simulated dataset, with mutational signatures names as colnames and channels as rownames.
#' @param true_exposures matrix of true exposures, for example obtained using simulations, with mutational signatures names as colnames and sample names as rownames. This is optional, and it useful to show which of the true_signature are actually present in the dataset.
#' @param outfile file name for plotting, please use .pdf file name extension. Can be omitted.
#' @return cosine similarity and signatures match between true and estimated signature names
#' @export
evaluatePerformanceSignatureSimilarity <- function(true_signatures,
                                                   estimated_signatures,
                                                   true_exposures = NULL,
                                                   outfile = NULL){
  # get the similarities
  distMatrix <- 1 - signature.tools.lib::computeCorrelationOfTwoSetsOfSigs(estimated_signatures,true_signatures)
  dummynames <- NULL

  if(ncol(estimated_signatures)!=ncol(true_signatures)){
    # message("[warning evaluateSignatureSimilarity] number of estimated signatures and truth signatures is not the same. Eliminating least similar signatures until same number is reached.")
    if (ncol(estimated_signatures)>ncol(true_signatures)){
      # need to remove estimated signatures
      ndiff <- ncol(estimated_signatures)-ncol(true_signatures)
      # add dummy variables
      dummynames <- (1:ndiff)+ncol(true_signatures)
      distMatrix <- cbind(distMatrix,matrix(data = 1,
                                            nrow = ncol(estimated_signatures),
                                            ncol = ndiff,
                                            dimnames = list(colnames(estimated_signatures),
                                                            dummynames)))

    }else if (ncol(estimated_signatures)<ncol(true_signatures)){
      # need to remove true_signatures signatures
      ndiff <- ncol(true_signatures)-ncol(estimated_signatures)
      # add dummy variables
      dummynames <- (1:ndiff)+ncol(estimated_signatures)
      distMatrix <- rbind(distMatrix,matrix(data = 1,
                                            nrow = ndiff,
                                            ncol=ncol(true_signatures),
                                            dimnames = list(dummynames,
                                                            colnames(true_signatures))))
    }
  }

  # find match
  res_match <- lpSolve::lp.assign(as.matrix(distMatrix))
  match_is <- apply(res_match$solution,1,which.max)
  # get cos sim
  matchCS <- 1 - apply(distMatrix*res_match$solution,1,max)

  unmatched <- matchCS < 0.85
  matchCS[unmatched] <- NA

  #return the match
  res <- list()
  res$matchTable <- data.frame(`estimated signatures`=rownames(distMatrix),
                               `matched true signatures`=colnames(distMatrix)[match_is],
                               `cosine similarity`=matchCS,stringsAsFactors = F,check.names = F)

  # still need to fix for dummy variables
  if (ncol(estimated_signatures)>ncol(true_signatures)){
    # remove dummy variables
    whereToRemove <- res$matchTable[,"matched true signatures"] %in% dummynames
    res$matchTable[whereToRemove,"matched true signatures"] <- NA
  }else if (ncol(estimated_signatures)<ncol(true_signatures)){
    # remove dummy variables
    whereToRemove <- res$matchTable[,"estimated signatures"] %in% dummynames
    res$matchTable[whereToRemove,"estimated signatures"] <- NA
  }

  res$minCosSim <- min(matchCS,na.rm = T)
  res$averageCosSim <- mean(matchCS,na.rm = T)

  # if possible annotate which signature are actually present in the dataset
  if(!is.null(true_exposures)){
    signames <- colnames(true_exposures)[apply(true_exposures,2,sum)>0]
    res$matchTable[,"true signature in dataset"] <- FALSE
    res$matchTable[res$matchTable$`matched true signatures` %in% signames,"true signature in dataset"] <- TRUE
  }

  # if needed, plot before returning
  if(!is.null(outfile)){
    plotPerformanceSignatures(matchTable = res$matchTable,
                              outfile = outfile)
  }

  return(res)
}

#' evaluate performance signature similarity list
#'
#' Compare a signatures matrix of true signatures with a list of matrices of estimated signatures, and
#' compute a match between the signatures and their cosine similarity.
#' For each matrix of estimated signatures, a match implies
#' an optimal assignment that maximises the cosine similarity while allowing each true_signature
#' to match at most one estimated signature and viceversa.
#' If the number of signatures differ, the least similar signatures are removed until the number
#' of signatures is the same
#' Results are plotted to a file if outfile is specified.
#'
#' @param true_signatures matrix of true signatures, for example obtained using simulations, with mutational signatures names as colnames and channels as rownames
#' @param estimated_signatures_list list of matrices of estimated signatures, for example obtained using a signature extraction on a simulated dataset, with mutational signatures names as colnames and channels as rownames. Use a named list, so that the names can be used in he data structure returned and plots.
#' @param true_exposures matrix of true exposures, for example obtained using simulations, with mutational signatures names as colnames and sample names as rownames. This is optional, and it useful to show which of the true_signature are actually present in the dataset.
#' @param outfile file name for plotting, please use .pdf file name extension. Can be omitted.
#' @return cosine similarity and signatures match between true and estimated signature names
#' @export
evaluatePerformanceSignatureSimilarityList <- function(true_signatures,
                                                       estimated_signatures_list,
                                                       true_exposures = NULL,
                                                       outfile = NULL){
  groupNames <- names(estimated_signatures_list)
  perfList <- list()
  for (g in groupNames) {
    perfList[[g]] <- evaluatePerformanceSignatureSimilarity(true_signatures,
                                                            estimated_signatures_list[[g]],
                                                            true_exposures = true_exposures)
  }

  signames <- c()
  for (g in groupNames){
    signames <- union(signames,perfList[[g]]$matchTable$`matched true signatures`)
  }
  # remove NA
  signames <- setdiff(signames,NA)
  returnTable <- data.frame(row.names = signames,check.names = F,stringsAsFactors = F)
  for (g in groupNames){
    tmpTable <- perfList[[g]]$matchTable
    returnTable[tmpTable$`matched true signatures`,g] <- tmpTable$`cosine similarity`
  }

  # if possible annotate which signature are actually present in the dataset
  if(!is.null(true_exposures)){
    sigsInData <- colnames(true_exposures)[apply(true_exposures,2,sum)>0]
    returnTable[,"true signature in dataset"] <- FALSE
    returnTable[sigsInData,"true signature in dataset"] <- TRUE
  }

  # if needed, plot before returning
  if(!is.null(outfile)){
    plotPerformanceSignaturesList(matchTable = returnTable,
                                  outfile = outfile)
  }

  returnObj <- list()
  returnObj$returnTable <- returnTable
  returnObj$perfList <- perfList

  return(returnObj)
}

plotPerformanceSignatures <- function(matchTable,
                                      outfile){
  dir.create(dirname(outfile),showWarnings = F,recursive = T)
  tmpMatrix <- matchTable[,3:ncol(matchTable),drop=F]
  row.names(tmpMatrix) <- matchTable$`matched true signatures`
  whichInDataset <- NULL
  if("true signature in dataset" %in% colnames(matchTable)){
    whichInDataset <- matchTable$`matched true signatures`[matchTable$`true signature in dataset`]
    tmpMatrix <- tmpMatrix[,-which(colnames(tmpMatrix)=="true signature in dataset"),drop=F]
  }

  sigsToPlot <- row.names(tmpMatrix)[apply(tmpMatrix,1,function(x) any(!is.na(x)))]
  if("true signature in dataset" %in% colnames(matchTable)){
    sigsToPlot <- union(sigsToPlot,whichInDataset)
  }
  tmpMatrix <- tmpMatrix[sigsToPlot,,drop=F]
  plotMatrix(tmpMatrix,
             output_file = outfile,
             thresholdMark = 0.95)
}


plotPerformanceSignaturesList <- function(matchTable,
                                          outfile){
  dir.create(dirname(outfile),showWarnings = F,recursive = T)
  tmpMatrix <- matchTable
  whichInDataset <- NULL
  if("true signature in dataset" %in% colnames(matchTable)){
    whichInDataset <- rownames(matchTable)[matchTable$`true signature in dataset`]
    tmpMatrix <- tmpMatrix[,-which(colnames(tmpMatrix)=="true signature in dataset"),drop=F]
  }

  sigsToPlot <- row.names(tmpMatrix)[apply(tmpMatrix,1,function(x) any(!is.na(x)))]
  if("true signature in dataset" %in% colnames(matchTable)){
    sigsToPlot <- union(sigsToPlot,whichInDataset)
  }
  tmpMatrix <- tmpMatrix[sigsToPlot,,drop=F]
  plotMatrix(tmpMatrix,
             output_file = outfile,
             thresholdMark = 0.95)
}

#' update signature names using a match table
#'
#' This function is part of the workflow of simulating a dataset, testing a
#' signature extraction and evaluating the signatures/exposures accuracy.
#' The function evaluatePerformanceSignatureSimilarity returns a table
#' matching true signatures with estimated signatures. The updateSigNamesWithMatchTable
#' function can be used with the match table to replace the names of the
#' estimated signatures with those of the true signatures, for example in the
#' matrix of estimated exposures. See example below.
#'
#' @param signames vector of signature names, typically estimated signatures, that need to be replaced
#' @param matchTable table with matched true and estimated signature names
#' @return updated signames vector
#' @export
#' @examples
#' sigsPerf <- evaluatePerformanceSignatureSimilarity(true_signatures,
#'                                                    estimated_signatures)
#' updatedNames <- updateSigNamesWithMatchTable(signames,
#'                                              sigsPerf$matchTable)
#'
updateSigNamesWithMatchTable <- function(signames,
                                         matchTable){
  tmpTable <- matchTable[complete.cases(matchTable),,drop=F]
  intnames <- intersect(signames,tmpTable$`estimated signatures`)
  for (i in intnames){
    signames[signames==i] <- tmpTable$`matched true signatures`[tmpTable$`estimated signatures`==i]
  }
  return(signames)
}
