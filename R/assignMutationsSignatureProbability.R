#' Assign Signature Probability to Mutations
#' 
#' Given an annotated set of mutations and an exposure matrix from signature fit, this function estimates the probability that each mutation generates from each signature.
#' 
#' @param sampleMutations list of annotated mutations obtained usually after building catalogues, for example using the functions vcfToSNVcatalogues or bedpeToRearrCatalogues
#' @param sampleSigsExposures matrix of exposures with only one row (one sample) and the exposures for each signature as columns
#' @param signatures mutational signatures matrix with the signatures used during signature fitting
#' @return matrix of mutations with an additional column containing the probabilities in a text format that can be expanded into a matrix using expandColumnToMatrix
#' @export
assignSignatureProbabilityToMutations <- function(sampleMutations,
                                                  sampleSigsExposures,
                                                  signatures,
                                                  verbose=TRUE){
  # some checks
  if(!nrow(sampleSigsExposures)==1){
    message("[error assignMutationsSignatureProbability] sampleSigsExposures should be a table with one row, ",
            "indicating signature exposures for one sample, however nrow=",nrow(sampleSigsExposures))
    return(sampleMutations)
  }
  # remove unassigned, if present
  if("unassigned" %in% colnames(sampleSigsExposures)){
    sampleSigsExposures <- sampleSigsExposures[,-which(colnames(sampleSigsExposures)=="unassigned"),drop=F]
  }
  # check that all exposures have a corresponding signature in the signatures table
  sigsAvailable <- colnames(sampleSigsExposures) %in% colnames(signatures)
  if(!all(sigsAvailable)){
    message("[error assignMutationsSignatureProbability] sampleSigsExposures contains exposures for signatures ",
            "that are not available in the provided signatures table. Signatures missing: ",
            paste(colnames(sampleSigsExposures)[!sigsAvailable],collapse = ", "),".")
    return(sampleMutations)
  }
  
  # attempt to identify the mutation type automatically
  muttype <- getTypeOfMutationsFromChannels(signatures)
  if(muttype=="subs"){
    channel_column <- "context"
  }else if(muttype=="DNV"){
    channel_column <- "channelAlexandrov"
  }else if(muttype=="rearr"){
    channel_column <- "catalogue.label"
  }else{
    channel_column <- "channel"
    if(verbose) message("[info assignMutationsSignatureProbability] signatures mutation type not recognised (are you using custom channels?) ",
                        "Attempting to use column \"channel\" in sampleMutations to match mutations with rownames of the signatures table. ",
                        "You can suppress this message using verbose=FALSE.")
  }
  
  # check that the required column is in the table and that the values match the signature
  if(!channel_column %in% colnames(sampleMutations)){
    message("[error assignMutationsSignatureProbability] sampleMutations table is missing the required \"",channel_column,"\" column.")
    return(sampleMutations)
  }
  # check that all channels map to signature channels
  mutsChannels <- unique(sampleMutations[,channel_column,drop=T])
  sigsChannels <- rownames(signatures)
  channelsAvailable <- mutsChannels %in% sigsChannels
  if(!all(channelsAvailable)){
    message("[error assignMutationsSignatureProbability] sampleMutations table contains mutation classes that are not among the rownames of the signatures table ",
            "The channels not found in the signature table are: ",paste(mutsChannels[!channelsAvailable],collapse = ", "),".")
    return(sampleMutations)
  }
  
  # OK now we are ready to calculate the probabilities
  # use only >0 exp and corresponding signatures
  sampleSigsExposures <- sampleSigsExposures[,sampleSigsExposures[1,]>0,drop=F]
  signatures <- signatures[,colnames(sampleSigsExposures),drop=F]
  # now if there are no positive exposures the result in compressed format will be an empty string
  if(ncol(sampleSigsExposures)==0){
    sampleMutations$sigsProb <- rep("",nrow(sampleMutations))
  }else{
    # scale sigs according to exposures
    weightedSigs <- matrix(as.matrix(sampleSigsExposures),
                           ncol = ncol(sampleSigsExposures),
                           nrow = nrow(signatures),
                           byrow = T) * as.matrix(signatures)
    # normalise by row to get channel probabilities
    channelProbs <-  as.matrix(weightedSigs/matrix(apply(weightedSigs,1,sum),
                                                   ncol = ncol(weightedSigs),
                                                   nrow = nrow(weightedSigs),
                                                   byrow = F))
    
    channelProbs[is.infinite(channelProbs)] <- 0
    channelProbs[is.nan(channelProbs)] <- 0
    sampleMutations$sigsProb <- rep("",nrow(sampleMutations))
    for(x in mutsChannels){
      probVect <- channelProbs[x,,drop=F]
      sampleMutations$sigsProb[sampleMutations[,channel_column]==x] <- paste(paste(colnames(probVect),collapse = ":"),paste(probVect[1,,drop=T],collapse = ":"),sep = ";")
    }
    
  }
  
  # return updated matrix
  return(sampleMutations)
}


#' Expand columns to matrix
#' 
#' @param dataTable table containing a column to expand. Each cell to expand will be a list of pairs, where each element of the list is separated by ";" and each pair is separated by ":". Pairs indicate colname and value. For example, if the value in a row is C1:V1;C2:V2;C4:V4, this will be expanded into columns C1, C2, C4 and values for that row V1, V2, V4.   
#' @param targetColname colname of the column to expand
#' @return dataTable with expanded column
#' @export
expandColumnToMatrix <- function(dataTable,
                                 targetColname){
  newcolnames <- c()
  for(i in 1:nrow(dataTable)){
    x <- dataTable[i,targetColname]
    if(x!="") newcolnames <- unique(c(newcolnames,strsplit(strsplit(x,split = ";")[[1]][1],split = ":")[[1]]))
  }
  
  newtable <- NULL
  if(length(newcolnames)>0){
    newtable <- as.data.frame(matrix(0,
                                     nrow = nrow(dataTable),
                                     ncol = length(newcolnames),
                                     dimnames = list(rownames(dataTable),newcolnames)),
                              stringsAsFactors = F)
    for(i in 1:nrow(dataTable)){
      x <- dataTable[i,targetColname]
      if(x!="") {
        sigsNames <- strsplit(strsplit(x,split = ";")[[1]][1],split = ":")[[1]]
        probVals <- as.numeric(strsplit(strsplit(x,split = ";")[[1]][2],split = ":")[[1]])
        newtable[i,sigsNames] <- probVals
      }
    }
    colnames(newtable) <- paste0("sigProb.",colnames(newtable))
  }
  dataTable <- dataTable[,-which(colnames(dataTable)==targetColname),drop=F]
  if(!is.null(newtable)) dataTable <- cbind(dataTable,newtable)
  return(dataTable)
}

