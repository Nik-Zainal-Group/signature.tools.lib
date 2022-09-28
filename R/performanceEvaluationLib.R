

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

#' @export
evaluatePerformanceExposures <- function(true_exposures,
                                         estimated_exposures,
                                         commonNames = NULL,
                                         rareNames = NULL,
                                         outfile = NULL){
  if(!is.null(commonNames)){
    resPerfCommon <- evaluatePerformanceExposuresCore(true_exposures[,commonNames,drop=F],
                                                      estimated_exposures[,intersect(commonNames,colnames(estimated_exposures)),drop=F])
  }else{
    intersectionNames <- intersect(colnames(true_exposures),colnames(estimated_exposures))
    resPerfCommon <- evaluatePerformanceExposuresCore(true_exposures[,intersectionNames,drop=F],
                                                      estimated_exposures[,intersectionNames,drop=F])
  }

  if(!is.null(rareNames)){
    resPerfRare <- evaluatePerformanceExposuresCore(true_exposures[,rareNames,drop=F],
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
  # perfColNames <- c()
  # perfTableNames <- c()
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
  nmetrics <- nrow(perfTables$resPerfCommon)
  if(is.null(perfTables$resPerfRare)) {
    plotRareSigs <- FALSE
  }else if (!(ncol(perfTables$resPerfRare)==0)){
    plotRareSigs <- TRUE
  }
  allcolours <- c("#F3C300","#875692","#F38400","#A1CAF1","#BE0032","#C2B280","#848482","#008856","#E68FAC","#0067A5","#F99379","#604E97","#F6A600","#B3446C","#DCD300")
  
  pointsize <- 16
  maxNamesLength <- max(sapply(colnames(perfTables$resPerfCommon),function(x) strwidth(x,units = "inch",ps = pointsize),USE.NAMES = F))
  
  plotHeight <- 5*nmetrics
  plotWidth <- ifelse(plotRareSigs,12,6)
  
  cairo_pdf(filename = outfile,height = plotHeight,width = plotWidth,pointsize = pointsize)
  par(mfrow=c(nmetrics,ifelse(plotRareSigs,2,1)),mai=c(maxNamesLength+0.1,1,1,1))
  for(i in 1:nmetrics){
    barplot(unlist(perfTables$resPerfCommon[i,,drop=T]),beside = T,
            names.arg = colnames(perfTables$resPerfCommon),
            col = allcolours[1:ncol(perfTables$resPerfCommon)],border = NA,
            main = rownames(perfTables$resPerfCommon)[i],las=2)
    if(plotRareSigs) barplot(unlist(perfTables$resPerfRare[i,,drop=T]),beside = T,
                             names.arg = colnames(perfTables$resPerfRare),
                             col = allcolours[1:ncol(perfTables$resPerfRare)],border = NA,
                             main = rownames(perfTables$resPerfRare)[i],las=2)
  }
  dev.off()
  return(perfTables)
}


#' @export
evaluatePerformanceSignatureSimilarity <- function(true_signatures,
                                                   estimated_signatures,
                                                   outfile = NULL){
  # get the similarities
  distMatrix <- 1 - signature.tools.lib::computeCorrelationOfTwoSetsOfSigs(estimated_signatures,true_signatures)
  
  if(ncol(estimated_signatures)!=ncol(true_signatures)){
    message("[warning evaluateSignatureSimilarity] number of estimated signatures and truth signatures is not the same. Eliminating least similar signatures until same number is reached.")
    if (ncol(estimated_signatures)>ncol(true_signatures)){
      # need to remove estimated signatures
      minsigdist <- apply(distMatrix,1,min)
      ndiff <- ncol(estimated_signatures)-ncol(true_signatures)
      for (i in 1:ndiff){
        ni <- which.max(minsigdist)
        minsigdist <- minsigdist[-ni]
        distMatrix <- distMatrix[-ni,]
      }
    }else if (ncol(estimated_signatures)<ncol(true_signatures)){
      # need to remove true_signatures signatures
      minsigdist <- apply(distMatrix,2,min)
      ndiff <- ncol(true_signatures)-ncol(estimated_signatures)
      for (i in 1:ndiff){
        ni <- which.max(minsigdist)
        minsigdist <- minsigdist[-ni]
        distMatrix <- distMatrix[,-ni]
      }
    }
  }
  
  # find match
  res_match <- lpSolve::lp.assign(as.matrix(distMatrix))
  match_is <- apply(res_match$solution,1,which.max)
  # get cos sim
  matchCS <- 1 - apply(distMatrix*res_match$solution,1,max)
  
  # check unmatched
  unmatched_estimated <- setdiff(colnames(estimated_signatures),rownames(distMatrix))
  unmatched_truth <- setdiff(colnames(true_signatures),colnames(distMatrix))
  
  #return the match
  res <- list()
  res$matchTable <- data.frame(`estimated signatures`=rownames(distMatrix),
                               `matched true signatures`=colnames(distMatrix)[match_is],
                               `cosine similarity`=matchCS,stringsAsFactors = F,check.names = F)
  # add other signatures unmatched
  if(length(unmatched_estimated)>0){
    res$matchTable <- rbind(res$matchTable,data.frame(`estimated signatures`=unmatched_estimated,
                                                      `matched true signatures`=rep(NA,length(unmatched_estimated)),
                                                      `cosine similarity`=rep(NA,length(unmatched_estimated)),
                                                      stringsAsFactors = F,check.names = F))
  }else if (length(unmatched_truth)>0){
    res$matchTable <- rbind(res$matchTable,data.frame(`estimated signatures`=rep(NA,length(unmatched_truth)),
                                                      `matched true signatures`=unmatched_truth,
                                                      `cosine similarity`=rep(NA,length(unmatched_truth)),
                                                      stringsAsFactors = F,check.names = F))
  }
  
  res$minCosSim <- min(matchCS)
  res$averageCosSim <- mean(matchCS)
  
  # if needed, plot before returning
  if(!is.null(outfile)){
    dir.create(dirname(outfile),showWarnings = F,recursive = T)
    tmpMatrix <- res$matchTable[,3:ncol(res$matchTable),drop=F]
    row.names(tmpMatrix) <- res$matchTable$`matched true signatures`
    plotMatrix(tmpMatrix,
               output_file = outfile,
               thresholdMark = 0.95)
  }
  
  return(res)
}

#' @export
evaluatePerformanceSignatureSimilarityList <- function(true_signatures,
                                                       estimated_signatures_list,
                                                       outfile = NULL){
  groupNames <- names(estimated_signatures_list)
  perfList <- list()
  for (g in groupNames) {
    perfList[[g]] <- evaluatePerformanceSignatureSimilarity(true_signatures,
                                                            estimated_signatures_list[[g]])
  }
  
  signames <- c()
  for (g in groupNames){
    signames <- union(signames,perfList[[g]]$matchTable$`matched true signatures`)
  }
  # remove NA
  signames <- setdiff(signames,NA)
  returnTable <- data.frame(row.names = signames,check.names = F,stringsAsFactors = F)
  for (g in groupNames){
    tmpTable <- perfList[[g]]$matchTable[complete.cases(perfList[[g]]$matchTable),,drop=F]
    returnTable[tmpTable$`matched true signatures`,g] <- tmpTable$`cosine similarity`
  }
  
  # if needed, plot before returning
  if(!is.null(outfile)){
    dir.create(dirname(outfile),showWarnings = F,recursive = T)
    plotMatrix(returnTable,
               output_file = outfile,
               thresholdMark = 0.95)
  }
  
  returnObj <- list()
  returnObj$returnTable <- returnTable
  returnObj$perfList <- perfList
  
  return(returnObj)
}


