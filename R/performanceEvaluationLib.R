

#' @export
evaluatePerformanceExposures <- function(true_exposures,
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
evaluatePerformanceExposuresCommonRare <- function(true_exposures,
                                                   estimated_exposures,
                                                   commonNames,
                                                   rareNames){
  resPerfCommon <- evaluatePerformanceExposures(true_exposures[,commonNames,drop=F],
                                                estimated_exposures[,intersect(commonNames,colnames(estimated_exposures)),drop=F])
  resPerfRare <- evaluatePerformanceExposures(true_exposures[,rareNames,drop=F],
                                                estimated_exposures[,intersect(rareNames,colnames(estimated_exposures)),drop=F])
  return(list(resPerfCommon=resPerfCommon,
              resPerfRare=resPerfRare))
}

#' @export
evaluatePerformanceSignatureSimilarity <- function(true_signatures,
                                                   estimated_signatures){
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
  return(res)
}

#' @export
combinePerformanceSignatureSimilarity <- function(tablesList){
  signames <- c()
  for (n in names(tablesList)){
    signames <- union(signames,tablesList[[n]]$`matched true signatures`)
  }
  returnTable <- data.frame(row.names = signames,check.names = F,stringsAsFactors = F)
  for (n in names(tablesList)){
    returnTable[tablesList[[n]]$`matched true signatures`,n] <- tablesList[[n]]$`cosine similarity`
  }
  return(returnTable)
}

#' @export
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
