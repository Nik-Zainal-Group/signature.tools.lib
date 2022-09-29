
#' simpleCatalogueSimulation
#'
#' This is a simple function that produces simulated mutational catalgoues and exposures, given a set of mutational signatures.
#' Each catalogue is simulated as follows. The total number of mutations in a catalogue is sampled uniformely in the log space,
#' and the results exponetiated, so that most samples have total mutations closer to min_mutations than max_mutations.
#' exposures for common signatures and rare signatures (if any) are then sampled uniformely, and the result scaled so that the
#' sum of the mutations corresponds to the number sampled initially.
#'
#' @param commonSignatures mutational signatures matrix with the signatures to be considered common
#' @param rareSignatures mutational signatures matrix with the signatures to be considered rare, can be omitted
#' @param ngenomes how many mutational catalogues should be simulated
#' @param ncommon_signatures_per_sample how many common signatures should be present in each sample
#' @param nsamplesWrareSigs how many samples should have rare signatures
#' @param maxRareSigsPerSample the maxiumum number of rare signature that can be present in a sample
#' @param min_mutations minimum number of mutations in a sample
#' @param max_mutations maximum number of mutations in a sample
#' @param randomSeed set a random seed for reproducibility
#' @return catalogues, signatures and exposures of simulated data
#' @export
simpleCatalogueSimulation <- function(commonSignatures,
                                      rareSignatures=NULL,
                                      ngenomes=50,
                                      ncommon_signatures_per_sample=5,
                                      nsamplesWrareSigs=10,
                                      maxRareSigsPerSample = 1,
                                      min_mutations=1000,
                                      max_mutations=50000,
                                      randomSeed = NULL){

  if(!is.null(randomSeed)){
    set.seed(randomSeed)
  }

  nmutations <- exp(runif(ngenomes,log(min_mutations),log(max_mutations)))
  nsignatures <- ncol(commonSignatures)
  exposures <- matrix(runif(ngenomes*nsignatures),nrow = nsignatures)

  #make more sparse
  if (nsignatures>ncommon_signatures_per_sample){
    for (i in 1:ncol(exposures)){
      selected_rows <- sample(1:nsignatures,nsignatures - ncommon_signatures_per_sample)
      exposures[selected_rows,i] <- 0
    }
  }

  allSignatures <- commonSignatures

  if(!is.null(rareSignatures)){
    allSignatures <- cbind(allSignatures,rareSignatures)
    nsignatures <- ncol(allSignatures)
    whichSamples <- sample(1:ngenomes,size = nsamplesWrareSigs,replace = F)
    rareExposures <- matrix(0,nrow = ncol(rareSignatures),ncol = ngenomes)
    for(i in 1:length(whichSamples)){
      howmanyraresigs <- sample(maxRareSigsPerSample,size = 1)
      whichSignatures <- sample(1:ncol(rareSignatures),size = howmanyraresigs,replace = T)
      rareExposures[whichSignatures,whichSamples[i]] <- runif(howmanyraresigs)
    }
    exposures <- rbind(exposures,rareExposures)
  }

  exposures <- exposures/matrix(rep(apply(exposures,2,sum),nsignatures),nrow = nsignatures,byrow = TRUE)
  exposures <- round(exposures*matrix(rep(nmutations,nsignatures),nrow = nsignatures,byrow = TRUE))

  #compose simulated catalogues
  catalogues <- mutationSampling(allSignatures,exposures)

  row.names(catalogues) <- row.names(allSignatures)
  colnames(catalogues) <- paste0("Sample_",1:ngenomes)

  colnames(exposures) <- paste0("Sample_",1:ngenomes)
  rownames(exposures) <- colnames(allSignatures)

  #add poisson noise
  for (i in nrow(catalogues)){
    catalogues[i,] <- rpois(ngenomes,catalogues[i,])
  }

  return(list(catalogues=catalogues,
              exposures=exposures,
              signatures=allSignatures))
}

#' mutationSampling
#'
#' @param allSignatures mutational signatures matrix with the signatures to sample
#' @param exposures exposures with how many muts for each signatures in each sample

#' @return sampled catalogues
#' @export
mutationSampling <- function(allSignatures,exposures){
  sampled <- matrix(0,nrow = nrow(allSignatures),ncol = ncol(exposures),dimnames = list(rownames(allSignatures),colnames(exposures)))
  for(j in 1:ncol(exposures)){
    for(i in 1:nrow(exposures)){
      if(exposures[i,j]>0){
        tmpSample <- table(sample(1:nrow(allSignatures),size = exposures[i,j],replace = T,prob = allSignatures[,i]))
        sampled[as.numeric(names(tmpSample)),j] <- sampled[as.numeric(names(tmpSample)),j] + tmpSample
      }
    }
  }
  return(sampled)
}


#' Save simulated data
#'
#' Function that can be used to save to file a dataset simulated with the
#' simpleCatalogueSimulation function.
#'
#' @param simObj object returned by the simpleCatalogueSimulation function
#' @param simDir directory where the simulated data should be stored
#' @param saveAsPlainText if TRUE data will be saved as plain text, in tab separated table format, while if FALSE data will be saved as R object.
#' @export
saveSimulatedData <- function(simObj,
                              simDir,
                              saveAsPlainText = T){
  if(!all(names(simObj) %in% c("catalogues","exposures","signatures"))){
    message("[error saveSimulatedData] simObj does not appear to be a valid simulated data object as returned from the simpleCatalogueSimulation function. Nothing saved.")
    return(NULL)
  }
  dir.create(simDir,showWarnings = F,recursive = T)
  if(!saveAsPlainText){
    save(file = paste0(simDir,"/simdata.rData"),simObj)
  }else{
    signature.tools.lib::writeTable(simObj$catalogues,paste0(simDir,"/catalogues.tsv"))
    signature.tools.lib::writeTable(simObj$exposures,paste0(simDir,"/exposures.tsv"))
    signature.tools.lib::writeTable(simObj$signatures,paste0(simDir,"/signatures.tsv"))
  }
}

#' Load simulated data
#'
#' Function that loads simulated data that was previously saved with the saveSimulatedData
#' function. The function can load  both data saved as R object and as plain text,
#' identifying the correct format automatically.
#'
#' @param simDir directory where the simulated data was saved
#' @return simulated data object
#' @export
loadSimulatedData <- function(simDir){
  simdatafile <- paste0(simDir,"/simdata.rData")
  simtxtfile <- paste0(simDir,"/catalogues.tsv")
  if(file.exists(simdatafile)){
    load(simdatafile)
  }else if(file.exists(simtxtfile)){
    simObj <- list()
    simObj$catalogues <- signature.tools.lib::readTable(paste0(simDir,"/catalogues.tsv"))
    simObj$exposures <- signature.tools.lib::readTable(paste0(simDir,"/exposures.tsv"))
    simObj$signatures <- signature.tools.lib::readTable(paste0(simDir,"/signatures.tsv"))
  }else{
    message("[error loadSimulatedData] simulated data not found in directory ",simDir)
    simObj <- NULL
  }

  return(simObj)
}

