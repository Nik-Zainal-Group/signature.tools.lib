
#' Mutational Signatures Fit (deprecated)
#' 
#' This function is deprecated. You can still use it, but we advise to use the function Fit instead, 
#' which provides a unified interface for basic signature fit with/without bootstrap.
#' Fit a given set of mutational signatures into mutational catalogues to extimate the activty/exposure of each of the given signatures in the catalogues.
#' 
#' @param cat catalogue matrix, patients as columns, channels as rows
#' @param signature_data_matrix signatures, signatures as columns, channels as rows
#' @param method KLD or NNLS or SA
#' @param bf_method bleeding filter method, one of KLD or CosSim, only if bleeding filter is used (alpha>-1)
#' @param alpha set alpha to -1 to avoid Bleeding Filter
#' @param doRound round the exposures to the closest integer
#' @param verbose use FALSE to suppress messages
#' @param n_sa_iter set max Simulated Annealing iterations if method==SA
#' @param showDeprecated set to FALSE to switch off the deprecated warning messsage
#' @return returns the activities/exposures of the signatures in the given sample
#' @keywords mutational signatures fit
#' @export
#' @examples
#' res.exp <- SignatureFit(catalogues,signature_data_matrix)
SignatureFit <- function(cat, #catalogue, patients as columns, channels as rows
                         signature_data_matrix,  #signatures, signatures as columns, channels as rows
                         method = "KLD", #KLD or NNLS or SA
                         bf_method = "CosSim", #KLD or CosSim
                         alpha = -1, #set alpha to -1 to avoid Bleeding Filter
                         doRound = TRUE, #round the exposures to the closest integer
                         verbose = TRUE, #use FALSE to suppress messages
                         n_sa_iter = 500,
                         showDeprecated = TRUE){ 
  if(showDeprecated) message("[warning SignatureFit] SignatureFit is deprecated, please use Fit instead with useBootstrap=FALSE. You can turn off this warning with showDeprecated=FALSE")
  
  if(method=="KLD"){
    if(verbose) message("[info SignatureFit] SignatureFit, objective function: KLD")
    
    #Fit using the given signatures
    nnlm_res <- NNLM::nnlm(as.matrix(signature_data_matrix),as.matrix(cat),loss = "mkl",method = "lee")
    fit_KLD <- KLD(cat,as.matrix(signature_data_matrix) %*% nnlm_res$coefficients)
    exposures <- nnlm_res$coefficients
    if(verbose) message("[info SignatureFit] Optimisation terminated, KLD=",fit_KLD)
    if(alpha >= 0){ #negative alpha skips Bleeding Filter
      #apply bleeding filter
      if(bf_method=="CosSim"){
        if(verbose) message("[info SignatureFit] Applying Bleeding Filter (Cosine Similarity) with alpha=",alpha*100,"%")
      }else if (bf_method=="KLD"){
        if(verbose) message("[info SignatureFit] Applying Bleeding Filter (KLD) with alpha=",alpha*100,"%")
      }
      for(i in 1:ncol(nnlm_res$coefficients)){
        if(bf_method=="CosSim"){
          exposures[,i] <- bleedingFilter(e = nnlm_res$coefficients[,i],
                                             sig = signature_data_matrix,
                                             sample = cat[,i],
                                             alpha = alpha)
        }else if (bf_method=="KLD"){
          exposures[,i] <- bleedingFilterKLD(e = nnlm_res$coefficients[,i],
                                             sig = signature_data_matrix,
                                             sample = cat[,i],
                                             alpha = alpha)
        }
      }
      fit_KLD_afterBF <- KLD(cat,as.matrix(signature_data_matrix) %*% exposures)
      if(verbose) message("[info SignatureFit] New fit after Bleeding Filter, KLD=",fit_KLD_afterBF," (",sprintf("%.3f",(fit_KLD_afterBF - fit_KLD)/fit_KLD*100),"% increase)")
    }
  }else if(method=="SA"){
    # library(GenSA)
    # library(foreach)
    # library(doParallel)
    # library(doMC)
    
    doParallel::registerDoParallel(4)
    
    if(verbose) message("[info SignatureFit] SignatureFit, objective function: SA")
    
    sig <- signature_data_matrix
    
    
    #---- copy from Sandro's code below
    exp_list <- foreach::foreach(j=1:ncol(cat)) %dopar% {	
      curr_cat <- as.numeric(cat[,j])
      nmut <- sum(curr_cat)
      exp_tmp <- rep(0, ncol(sig))
      names(exp_tmp) <- names(sig)
      if(nmut>0){
        if(verbose) message("[info SignatureFit] Analyzing " ,  j, " of ", ncol(cat), " sample name ", colnames(cat[j]))
        
        ## Compute the start solution for the sim. annealing
        ss <- startSolution(sig, curr_cat)
        
        ## Compute the exposure by using the sim. annealing
        out <- GenSA::GenSA(par=as.numeric(ss), lower=rep(0, ncol(sig)), upper=rep(nmut, ncol(sig)), fn=objSimAnnelaingFunction, control=list(maxit=n_sa_iter), xsample=curr_cat, xsignature=sig)
        
        exp_tmp <- as.numeric(out$par)
        if(sum(exp_tmp)==0){
          exp_tmp <- ss
        }else{
          exp_tmp <- (exp_tmp/sum(exp_tmp))*nmut
        }
        names(exp_tmp) <- names(ss)
        
        ## Apply the bleeding Filter to remove the 'unnecessary' signatures
        if(bf_method=="CosSim"){
          if(verbose) message("[info SignatureFit] Applying Bleeding Filter (Cosine Similarity) with alpha=",alpha*100,"%")
          exp_tmp <-  bleedingFilter(exp_tmp, sig, curr_cat, alpha)
        }else if (bf_method=="KLD"){
          if(verbose) message("[info SignatureFit] Applying Bleeding Filter (KLD) with alpha=",alpha*100,"%")
          exp_tmp <-  bleedingFilterKLD(exp_tmp, sig, curr_cat, alpha)
        }
        
      }
      exp_tmp
    }
    
    exposures <- matrix(unlist(exp_list), ncol=ncol(sig), byrow=T)
    rownames(exposures) <- colnames(cat)
    colnames(exposures) <- colnames(sig)
    exposures <- t(exposures)
    #----- end of Sandro's code
    

    
  }else if(method=="NNLS"){
    if(verbose) message("[info SignatureFit] SignatureFit, method: NNLS")
    
    exposures <- matrix(NA,ncol = ncol(cat),nrow = ncol(signature_data_matrix))
    colnames(exposures) <- colnames(cat)
    rownames(exposures) <- colnames(signature_data_matrix)
    for (i in 1:ncol(cat)){
      q <- as.vector(cat[,i])
      exp_NNLS <- nnls::nnls(as.matrix(signature_data_matrix),q)
      exposures[,i] <- exp_NNLS$x
    }
    
    fit_KLD <- KLD(cat,as.matrix(signature_data_matrix) %*% exposures)
    if(verbose) message("[info SignatureFit] Optimisation terminated, KLD=",fit_KLD)
    
    if(alpha >= 0){ #negative alpha skips Bleeding Filter
      #apply bleeding filter
      if(bf_method=="CosSim"){
        if(verbose) message("[info SignatureFit] Applying Bleeding Filter (Cosine Similarity) with alpha=",alpha*100,"%")
      }else if (bf_method=="KLD"){
        if(verbose) message("[info SignatureFit] Applying Bleeding Filter (KLD) with alpha=",alpha*100,"%")
      }
      for (i in 1:ncol(exposures)){
        if(bf_method=="CosSim"){
          exposures[,i] <- bleedingFilter(e = exposures[,i],
                                          sig = signature_data_matrix,
                                          sample = cat[,i],
                                          alpha = alpha)
        }else if (bf_method=="KLD"){
          exposures[,i] <- bleedingFilterKLD(e = exposures[,i],
                                             sig = signature_data_matrix,
                                             sample = cat[,i],
                                             alpha = alpha)
        }
      }
    }
    fit_KLD_afterBF <- KLD(cat,as.matrix(signature_data_matrix) %*% exposures)
    if(verbose) message("[info SignatureFit] New fit after Bleeding Filter, KLD=",fit_KLD_afterBF," (",sprintf("%.3f",(fit_KLD_afterBF - fit_KLD)/fit_KLD*100),"% increase)")
  }else{
    stop("[error SignatureFit] Unknown method specified.")
  }
  exposures[exposures < .Machine$double.eps] <- 0
  if(doRound) exposures <- round(exposures)
  return(exposures)
}



#' Mutational Signatures Fit with Bootstrap (deprecated)
#' 
#' This function is deprecated. You can still use it, but we advise to use the function Fit instead, 
#' which provides a unified interface for basic signature fit with/without bootstrap.
#' Fit a given set of mutational signatures into mutational catalogues to extimate the 
#' activty/exposure of each of the given signatures in the catalogues. Implementation 
#' of method similar to Huang et al. 2017, Detecting presence of mutational signatures with 
#' confidence, which uses a bootstrap apporach to calculate the empirical probability 
#' of an exposure to be larger or equal to a given threshold (i.e. 5% of mutations of a sample). 
#' This probability can be used to decide which exposures to remove from the initial fit, 
#' thus increasing the sparsity of the exposures.
#' 
#' @param cat catalogue matrix, patients as columns, channels as rows
#' @param signature_data_matrix signatures, signatures as columns, channels as rows
#' @param nboot number of bootstraps to use, more bootstraps more accurate results
#' @param exposureFilterType use either fixedThreshold or giniScaledThreshold. When using fixedThreshold, exposures will be removed based on a fixed percentage with respect to the total number of mutations (threshold_percent will be used). When using giniScaledThreshold each signature will used a different threshold calculated as (1-Gini(signature))*giniThresholdScaling
#' @param threshold_percent threshold in percentage of total mutations in a sample, only exposures larger than threshold are considered
#' @param giniThresholdScaling scaling factor for the threshold type giniScaledThreshold, which is based on the Gini score of a signature
#' @param threshold_p.value p-value to determine whether an exposure is above the threshold_percent. In other words, this is the empirical probability that the exposure is lower than the threshold
#' @param method KLD or NNLS or SA
#' @param bf_method bleeding filter method, one of KLD or CosSim, only if bleeding filter is used (alpha>-1)
#' @param alpha set alpha to -1 to avoid Bleeding Filter
#' @param doRound round the exposures to the closest integer
#' @param nparallel to use parallel specify >1
#' @param verbose use FALSE to suppress messages
#' @param n_sa_iter set max Simulated Annealing iterations if method==SA
#' @param randomSeed set an integer random seed
#' @param showDeprecated set to FALSE to switch off the deprecated warning messsage
#' @return returns the activities/exposures of the signatures in the given sample and other information, such as p-values and exposures of individual bootstrap runs.
#' @keywords mutational signatures fit
#' @references Huang, X., Wojtowicz, D., & Przytycka, T. M. (2017). Detecting Presence Of Mutational Signatures In Cancer With Confidence. bioRxiv, (October). https://doi.org/10.1101/132597
#' @export
#' @examples
#' res <- SignatureFit_withBootstrap(catalogues,signature_data_matrix)
SignatureFit_withBootstrap <- function(cat, #catalogue, patients as columns, channels as rows
                                       signature_data_matrix, #signatures, signatures as columns, channels as rows
                                       nboot = 100, #number of bootstraps to use, more bootstraps more accurate results
                                       exposureFilterType = "fixedThreshold", # or "giniScaledThreshold"
                                       giniThresholdScaling = 10,
                                       threshold_percent = 5, #threshold in percentage of total mutations in a sample, only exposures larger than threshold are considered
                                       threshold_p.value = 0.05, #p-value to determine whether an exposure is above the threshold_percent. In other words, this is the empirical probability that the exposure is lower than the threshold
                                       method = "KLD", #KLD or SA, just don't use SA or you will wait forever, expecially with many bootstraps. SA is ~1000 times slower than KLD or NNLS
                                       bf_method = "CosSim", #KLD or CosSim, only used if alpha != -1
                                       alpha = -1, #set alpha to -1 to avoid Bleeding Filter
                                       verbose=TRUE, #use FALSE to suppress messages
                                       doRound = FALSE, #round the exposures to the closest integer
                                       nparallel=1, #to use parallel specify >1
                                       n_sa_iter = 500, #only used if  method = "SA"
                                       randomSeed = NULL,
                                       showDeprecated = TRUE){ 
  if(showDeprecated) message("[warning SignatureFit_withBootstrap] SignatureFit_withBootstrap is deprecated, please use Fit instead with useBootstrap=TRUE. You can turn off this warning with showDeprecated=FALSE")
  if(!is.null(randomSeed)){
    set.seed(randomSeed)
  }
  
  if (nparallel > 1){
    # library(foreach)
    # library(doParallel)
    # library(doMC)
    doParallel::registerDoParallel(nparallel)
    if(!is.null(randomSeed)){
      doRNG::registerDoRNG(randomSeed)
    }
    boot_list <- foreach::foreach(j=1:nboot) %dopar% {
      bootcat <- generateRandMuts(cat)
      SignatureFit(bootcat,signature_data_matrix,method,bf_method,alpha,verbose=verbose,doRound = doRound,n_sa_iter=n_sa_iter,showDeprecated = F)
    }
  }else{
    boot_list <- list()
    for(i in 1:nboot){
      bootcat <- generateRandMuts(cat)
      boot_list[[i]] <- SignatureFit(bootcat,signature_data_matrix,method,bf_method,alpha,verbose=verbose,doRound = doRound,n_sa_iter=n_sa_iter,showDeprecated = F)
    }
  }
  
  samples_list <- list()
  for(i in 1:ncol(cat)) {
    samples_list[[colnames(cat)[i]]] <- matrix(NA,ncol = nboot,nrow = ncol(signature_data_matrix))
    colnames(samples_list[[colnames(cat)[i]]]) <- 1:nboot
    row.names(samples_list[[colnames(cat)[i]]]) <- colnames(signature_data_matrix)
  }
  for(i in 1:nboot){
    for(j in 1:ncol(cat)){
      samples_list[[colnames(cat)[j]]][,i] <- boot_list[[i]][,j]
    }
  }
  E_median_notfiltered <- matrix(NA,nrow = ncol(signature_data_matrix),ncol = ncol(cat))
  E_median_filtered <- matrix(NA,nrow = ncol(signature_data_matrix),ncol = ncol(cat))
  E_p.values <- matrix(NA,nrow = ncol(signature_data_matrix),ncol = ncol(cat))
  colnames(E_median_notfiltered) <- colnames(cat)
  row.names(E_median_notfiltered) <- colnames(signature_data_matrix)
  colnames(E_median_filtered) <- colnames(cat)
  row.names(E_median_filtered) <- colnames(signature_data_matrix)
  colnames(E_p.values) <- colnames(cat)
  row.names(E_p.values) <- colnames(signature_data_matrix)
  #KLD error vector
  KLD_samples <- c()
  
  for(i in 1:ncol(cat)) {
    if(sum(cat[,i])>0){
      boots_perc <- samples_list[[colnames(cat)[i]]]/matrix(apply(samples_list[[colnames(cat)[i]]],2,sum),byrow = TRUE,nrow = nrow(samples_list[[colnames(cat)[i]]]),ncol = ncol(samples_list[[colnames(cat)[i]]]))*100
      
      
      if(exposureFilterType=="giniScaledThreshold"){
        sigInvGini <- 1 - apply(signature_data_matrix,2,giniCoeff)
        giniThresholdPerc <- giniThresholdScaling*sigInvGini
        # set to zero differently for each signature
        p.values <- c()
        for(j in 1:length(giniThresholdPerc)) p.values <- c(p.values,sum(boots_perc[j,]<=giniThresholdPerc[j])/nboot)
        names(p.values) <- colnames(signature_data_matrix)
      }else if(exposureFilterType=="fixedThreshold"){
        p.values <- apply(boots_perc <= threshold_percent,1,sum)/nboot
      }
      
      
      median_mut <- apply(samples_list[[colnames(cat)[i]]],1,median)
      E_median_notfiltered[,i] <- median_mut
      E_p.values[,i] <- p.values
      
      median_mut_perc <- median_mut/sum(median_mut)*100
      # plot(median_mut_perc)
      # abline(h=5)
      median_mut_perc[p.values > threshold_p.value] <- 0
      #below rescaling, not sure whether to use it or not. If not I have something like a residual
      # median_mut_perc <- median_mut_perc/sum(median_mut_perc)*100
      median_mut <- median_mut_perc/100*sum(cat[,i])
      # boxplot(t(samples_list[[1]]))
      # points(1:10,E[,1],col="red")
      # points(1:10,median_mut,col="green")
      E_median_filtered[,i] <- median_mut
      KLD_samples <- c(KLD_samples,KLD(cat[,i,drop=FALSE],as.matrix(signature_data_matrix) %*% E_median_filtered[,i,drop=FALSE]))
    }else{
      E_median_notfiltered[,i] <- 0
      E_p.values[,i] <- NA
      E_median_filtered[,i] <- 0
      KLD_samples <- c(KLD_samples,0)
    }
  }
  names(KLD_samples) <- colnames(cat)
  
  #unassigned mutations
  reconstructed_with_median <- round(as.matrix(signature_data_matrix) %*% as.matrix(E_median_filtered))
  unassigned_muts <- sapply(1:ncol(reconstructed_with_median),function(i) sum(cat[,i,drop=FALSE]) - sum(reconstructed_with_median[,i,drop=FALSE]))
  names(unassigned_muts) <- colnames(cat)
  
  res <- list()
  res$signature_data_matrix <- signature_data_matrix
  res$cat <- cat
  res$E_median_filtered <- E_median_filtered
  res$E_p.values <- E_p.values
  res$samples_list <- samples_list
  res$boot_list <- boot_list
  res$KLD_samples <- KLD_samples
  #need metadata
  res$threshold_percent <- threshold_percent
  res$threshold_p.value <- threshold_p.value
  res$nboots <- nboot
  res$method <- method
  res$unassigned_muts <- unassigned_muts
  res$alpha <- alpha
  res$bf_method <- bf_method
  res$n_sa_iter <- n_sa_iter
  return(res)
}

#' Mutational Signatures Fit with Bootstrap Analysis
#' 
#' This function is a wrapper for the function SignatureFit_withBootstrap_Analysis, which
#' produces several plots for each sample in the catalogues cat.
#' Fit a given set of mutational signatures into mutational catalogues to extimate the 
#' activty/exposure of each of the given signatures in the catalogues. Implementation 
#' of method similar to Huang 2017, Detecting presence of mutational signatures with 
#' confidence, which uses a bootstrap apporach to calculate the empirical probability 
#' of an exposure to be larger or equal to a given threshold (i.e. 5% of mutations of a sample). 
#' This probability can be used to decide which exposures to remove from the initial fit, 
#' thus increasing the sparsity of the exposures.
#' Note that SignatureFit_withBootstrap_Analysis will save the results of SignatureFit_withBootstrap
#' in the outdir directory using the R save() function. If SignatureFit_withBootstrap_Analysis
#' is rerun with the same setting, the saved file will be loaded to avoid rerunning the Signature Fit
#' and figures will be replotted.
#' 
#' @param outdir output directory for the analysis, remember to add '/' at the end
#' @param cat catalogue matrix, patients as columns, channels as rows
#' @param signature_data_matrix signatures, signatures as columns, channels as rows
#' @param nboot number of bootstraps to use, more bootstraps more accurate results
#' @param type_of_mutations either "subs", "rearr" or "generic" 
#' @param threshold_percent threshold in percentage of total mutations in a sample, only exposures larger than threshold are considered
#' @param threshold_p.value p-value to determine whether an exposure is above the threshold_percent. In other words, this is the empirical probability that the exposure is lower than the threshold
#' @param method KLD or NNLS or SA
#' @param bf_method bleeding filter method, one of KLD or CosSim, only if bleeding filter is used (alpha>-1)
#' @param alpha set alpha to -1 to avoid Bleeding Filter
#' @param doRound round the exposures to the closest integer
#' @param nparallel to use parallel specify >1
#' @param verbose use FALSE to suppress messages
#' @param n_sa_iter set max Simulated Annealing iterations if method==SA
#' @return returns the activities/exposures of the signatures in the given sample and other information, such as p-values and exposures of individual bootstrap runs.
#' @keywords mutational signatures fit
#' @export
#' @examples
#' res <- SignatureFit_withBootstrap_Analysis(catalogues,signature_data_matrix)
SignatureFit_withBootstrap_Analysis <- function(outdir, #output directory for the analysis, remember to add '/' at the end
                                                cat, #catalogue, patients as columns, channels as rows
                                       signature_data_matrix, #signatures, signatures as columns, channels as rows
                                       nboot = 100, #number of bootstraps to use, more bootstraps more accurate results
                                       type_of_mutations="subs", #use one of c("subs","rearr","generic","dnv")
                                       threshold_percent = 5, #threshold in percentage of total mutations in a sample, only exposures larger than threshold are considered
                                       threshold_p.value = 0.05, #p-value to determine whether an exposure is above the threshold_percent. In other words, this is the empirical probability that the exposure is lower than the threshold
                                       method = "KLD", #KLD or SA, just don't use SA or you will wait forever, expecially with many bootstraps. SA is ~1000 times slower than KLD or NNLS
                                       bf_method = "CosSim", #KLD or CosSim, only used if alpha != -1
                                       alpha = -1, #set alpha to -1 to avoid Bleeding Filter
                                       doRound = TRUE, #round the exposures to the closest integer
                                       nparallel=1, #to use parallel specify >1
                                       n_sa_iter = 500){  #only used if  method = "SA"
  # outdir <- "../results/sigfitbootstraptests/"
  dir.create(outdir,recursive = TRUE,showWarnings = FALSE)
  
  #begin by computing the sigfit bootstrap
  file_store <- paste0(outdir,"SigFit_withBootstrap_Summary_m",method,"_bfm",bf_method,"_alpha",alpha,"_tr",threshold_percent,"_p",threshold_p.value,".rData")
  if(file.exists(file_store)){
    load(file_store)
    message("[info SignatureFit_withBootstrap_Analysis] Bootstrap Signature Fits loaded from file")
  }else{
    res <- SignatureFit_withBootstrap(cat = cat,
                                          signature_data_matrix = signature_data_matrix,
                                          nboot = nboot,
                                          threshold_percent = threshold_percent,
                                          threshold_p.value = threshold_p.value,
                                          method = method,
                                          bf_method = bf_method,
                                          alpha = alpha,
                                          doRound = doRound,
                                          nparallel = nparallel,
                                          n_sa_iter = n_sa_iter,
                                          verbose = FALSE)
    save(file = file_store,res,nboot)
  }
  
  plot_SignatureFit_withBootstrap(outdir,res,type_of_mutations)
  
  return(res)
}

## This Function removes the 'unnecesary' signatures
#optimise w.r.t. cosine similarity
#alpha is the max cosine similarity that can be lost 
bleedingFilter <- function(e, sig,  sample, alpha){
  ## Compute the cosine similarity between the current solution and the catalogue
  #sim_smpl <- computeSimSample(e, sample, sig)
  sim_smpl <- as.matrix(sig) %*% e
  val <- cos_sim(sample, sim_smpl)
  
  e_orig <- e
  delta <- 0
  
  
  while(delta<=alpha && length(which(e>0))>1){
    
    pos <- which(e>0)
    e <- e[pos]
    sig <- sig[,pos]
    sim_m <- matrix(0, ncol(sig), ncol(sig))
    colnames(sim_m) <- names(e)
    rownames(sim_m) <- names(e)
    
    ## Move mutations across each pair of signatures and estimate the cosine similarity of the new solution		
    for(i in 1:ncol(sig)){
      for(j in 1:ncol(sig)){
        if(i!=j){
          e2 <- e
          e2[j] <- e2[j]+e2[i]
          e2[i] <- 0
          #sim_smpl <- computeSimSample(e2, sample, sig)
          sim_smpl <- as.matrix(sig) %*% e2
          sim_m[i,j] <- cos_sim(sample, sim_smpl)	
        }
      }
    }
    
    ## Extract the minimum delta
    delta <- val-max(sim_m)
    
    # If delta <= alpha accept the new solution
    if(delta<=alpha){
      pos <-  which(sim_m==max(sim_m), arr.ind=T)
      e[pos[1,2]] <- e[pos[1,2]]+e[pos[1,1]]
      e[pos[1,1]] <- 0
    }
  }
  
  e_orig[] <- 0
  e_orig[names(e)] <- e
  
  return(e_orig)
}

## This Function removes the 'unnecesary' signatures
#optimise w.r.t. KLD
#alpha is the max ratio of KLD that can be lost (e.g. 0.01 is 1% of original KLD)
bleedingFilterKLD <- function(e, sig,  sample, alpha){
  ## Compute the cosine similarity between the current solution and the catalogue
  #sim_smpl <- computeSimSample(e, sample, sig)
  sim_smpl <- as.matrix(sig) %*% e
  val <- KLD(sample, sim_smpl)
  alpha <- alpha*val
  
  e_orig <- e
  delta <- 0
  
  
  while(delta<=alpha && length(which(e>0))>1){
    
    pos <- which(e>0)
    e <- e[pos]
    sig <- sig[,pos]
    sim_m <- matrix(0, ncol(sig), ncol(sig))
    colnames(sim_m) <- names(e)
    rownames(sim_m) <- names(e)
    
    ## Move mutations across each pair of signatures and estimate the cosine similarity of the new solution		
    for(i in 1:ncol(sig)){
      for(j in 1:ncol(sig)){
        if(i!=j){
          e2 <- e
          e2[j] <- e2[j]+e2[i]
          e2[i] <- 0
          #sim_smpl <- computeSimSample(e2, sample, sig)
          sim_smpl <- as.matrix(sig) %*% e2
          sim_m[i,j] <- KLD(sample, sim_smpl)	
        }
      }
    }
    sim_m <- sim_m + diag(nrow(sim_m))*1e6
    ## Extract the minimum delta
    delta <- min(sim_m)-val
    
    # If delta <= alpha accept the new solution
    if(delta<=alpha){
      pos <-  which(sim_m==min(sim_m), arr.ind=T)
      e[pos[1,2]] <- e[pos[1,2]]+e[pos[1,1]]
      e[pos[1,1]] <- 0
    }
  }
  
  e_orig[] <- 0
  e_orig[names(e)] <- e
  
  return(e_orig)
}


## Function to generate the start solution for the sim. annelaing
startSolution <- function(sig, cat){
  summ <- apply(sig, 1, sum)
  out <- (sig/summ)*cat
  pos <- which(is.nan(as.matrix(out)), arr.ind=T)
  if(length(pos)>0){
    out[pos] <- 0
  }
  return(apply(out, 2, sum))
}

## The Objective Function fot the simulated annelaing
objSimAnnelaingFunction <- function(x, xsample, xsignature){
  sim_smpl <- rep(0, length(xsample))
  for(i in 1:ncol(xsignature)){
    sim_smpl  <- sim_smpl+(xsignature[,i]*x[i])
  }
  sum(abs(xsample-sim_smpl))
}

## Given the exposure and the probability matrix compute the Catalogue 
# computeSimSample <- function(x, xsample, xsignature){
#   sim_smpl <- rep(0, length(xsample))
#   for(i in 1:ncol(xsignature)){
#     sim_smpl  <- sim_smpl+(xsignature[,i]*x[i])
#   }
#   return(sim_smpl)
# }




#' @export
plotExposures <- function(exposures,
                          output_file=NULL){
  
  plotMatrix(dataMatrix = exposures,
             output_file = output_file,
             ndigitsafterzero = 0)
}

#' Export Signature Fit with bootstrap to JSON
#' 
#' Given a res file obtained from the SignatureFit_withBootstrap or 
#' SignatureFit_withBootstrap_Analysis function, export it to a set of 
#' JSON files that can be used for web visualisation
#' 
#' @param outdir output directory where the output JSON files should be saved, remember to put "/" at the end of the folder name.
#' @param res R object obtained from the SignatureFit_withBootstrap or SignatureFit_withBootstrap_Analysis function
#' @keywords mutational signatures fit JSON
#' @export
#' @examples
#' res <- SignatureFit_withBootstrap(cat = rnd_matrix,
#'                   signature_data_matrix = cosmic30,
#'                   nboot = 5,
#'                   threshold_percent = 0.1,
#'                   threshold_p.value = 0.1)
#' export_SignatureFit_withBootstrap_to_JSON("JSON_out/",res)
export_SignatureFit_withBootstrap_to_JSON <- function(outdir,res){
  
  dir.create(outdir,showWarnings = FALSE,recursive = TRUE)
  #plot consensus exposures file with metadata
  consensus_file <- paste0(outdir,"consensus.json")
  sink(file = consensus_file)
  cat("{\n")
  cat("\t\"nboots\": ",res$nboots,",\n",sep = "")
  cat("\t\"threshold_percent\": ",res$threshold_percent,",\n",sep = "")
  cat("\t\"threshold_p.value\": ",res$threshold_p.value,",\n",sep = "")
  cat("\t\"method\": \"",res$method,"\",\n",sep = "")
  cat("\t\"consensus\": {\n",sep = "")
  
  for(i in 1:ncol(res$E_median_filtered)){
    sname <- colnames(res$E_median_filtered)[i]
    cat("\t\t\"",sname,"\": {\n",sep = "")
    
    for(j in 1:nrow(res$E_median_filtered)){
      rname <- rownames(res$E_median_filtered)[j]
      cat("\t\t\t\"",rname,"\": ",res$E_median_filtered[j,i],sep = "")
      if(j<nrow(res$E_median_filtered)){
        cat(",\n")
      }else{
        cat("\n")
      }
    }
    
    cat("\t\t}")
    if(i<ncol(res$E_median_filtered)){
      cat(",\n")
    }else{
      cat("\n")
    }
  }
  
  cat("\t}\n")
  cat("}\n")
  sink()
  
  #plot bootstraps exposures file 
  boot_file <- paste0(outdir,"bootstraps.json")
  sink(file = boot_file)
  cat("[\n")
  
  for (b in 1:length(res$boot_list)){
    data_mat <- res$boot_list[[b]]
    
    cat("\t{\n")
    
    for(i in 1:ncol(data_mat)){
      sname <- colnames(data_mat)[i]
      cat("\t\t\"",sname,"\": {\n",sep = "")
      
      for(j in 1:nrow(data_mat)){
        rname <- rownames(data_mat)[j]
        cat("\t\t\t\"",rname,"\": ",data_mat[j,i],sep = "")
        if(j<nrow(data_mat)){
          cat(",\n")
        }else{
          cat("\n")
        }
      }
      
      cat("\t\t}")
      if(i<ncol(data_mat)){
        cat(",\n")
      }else{
        cat("\n")
      }
    }
    
    cat("\t}")
    if(b<length(res$boot_list)){
      cat(",\n")
    }else{
      cat("\n")
    }
    
  }
  cat("]\n")
  sink()
  
  #plot correlation of exposures file for each sample
  for (s in 1:ncol(res$E_median_filtered)){
    if (nrow(res$samples_list[[s]])>1){
      sname <- colnames(res$E_median_filtered)[s]
      data_mat <- cor(t(res$samples_list[[s]]),method = "spearman")
      data_mat[is.na(data_mat)] <- 0
      data_mat[row(data_mat)+(ncol(data_mat)-col(data_mat))>=ncol(data_mat)] <- 0
      
      cor_file <- paste0(outdir,sname,"_correlation.tsv")
      sink(file = cor_file)
      cat("Signature")
      for (j in colnames(data_mat)) cat("\t",j,sep = "")
      cat("\n")
      
      for(i in 1:nrow(data_mat)){
        cat(row.names(data_mat)[i])
        for (j in 1:ncol(data_mat)) cat("\t",data_mat[i,j],sep = "")
        cat("\n")
      }
      
      sink()
    }
  }

}

#' Distribution of Signatures in Samples
#' 
#' Given a catalogue of samples and an exposures table, compute the relative amount of each signature
#' in each sample and the unassigned mutations. Also cluster the samples with hierarchical clustering
#' with average linkage and order the samples according to the clustering. Optionally, plot to file.
#' 
#' @param fileout if specified, generate a plot, otherwise no plot is generated
#' @param catalogue original catalogue, channels as rows and samples as columns
#' @param exposures exposures/activities of signatures in each sample. Signatures as rows, samples as columns
#' @keywords unexplained samples
#' @return list with to objects: the matrix of the distribution of the signatures in the samples and the hierachical clustering object
#' @export
#' @examples
#' res <- SignatureFit_withBootstrap(cat = catalogue,
#'                   signature_data_matrix = cosmic30,
#'                   nboot = 5,
#'                   threshold_percent = 0.1,
#'                   threshold_p.value = 0.1)
#' distribution_object <- exposureDistributionBarplot(catalogue=catalogue,
#'                   exposures=res$E_median_filtered)
exposureDistributionBarplot <- function(fileout=NULL,catalogue,exposures){
  total_mut_in_exp <- apply(exposures,2,sum)
  total_mut_in_cat <- apply(catalogue,2,sum)
  unassigned_mut <- total_mut_in_cat - total_mut_in_exp
  unassigned_mut[unassigned_mut < 0] <- 0
  plot_matrix <- rbind(exposures,unassigned_mut)
  plot_matrix <- plot_matrix/matrix(rep(total_mut_in_cat,nrow(plot_matrix)),byrow = TRUE,nrow = nrow(plot_matrix))*100
  #cluster samples
  d <- dist(t(plot_matrix), method = "euclidean") # distance matrix
  fit <- hclust(d, method="average") 
  #order according to clustering
  plot_matrix <- plot_matrix[,fit$order]
  kelly_colors <- c('F2F3F4', '222222', 'F3C300', '875692', 'F38400', 'A1CAF1', 'BE0032', 
                    'C2B280', '848482', '008856', 'E68FAC', '0067A5', 'F99379', '604E97', 
                    'F6A600', 'B3446C', 'DCD300', '882D17', '8DB600', '654522', 'E25822', '2B3D26','CCCCCC','CCCCCC','CCCCCC')
  kelly_colors <- paste0("#",kelly_colors)
  if(!is.null(fileout)){
    jpeg(filename = fileout,width = max(1800,200+ncol(exposures)*3),height = 1000,res = 200)
    par(mar=c(4,4,3,5),mgp=c(1.5,0.5,0))
    barplot(plot_matrix,col = c(kelly_colors[1:nrow(exposures)],"grey"),
            border = NA,ylab = "% of mutation",xlab = "samples",xaxt="n",space=0)
    legend(x="topright",title = "Signatures",legend = c(rownames(exposures),"other"),fill = c(kelly_colors[1:nrow(exposures)],"grey"),bty = "n",inset = c(-0.08,0),xpd = TRUE)
    dev.off()
  }
  res <- list()
  res$distribution_matrix <- plot_matrix
  res$clustering_object <- fit
  return(res)
}


#' Estimate samples not fully explained by signature fit
#' 
#' Given a catalogue of samples, signatures and exposures, compute the sum of the absolute deviations (SAD)
#' between the original catalogue and the reconstructed samples (i.e. signatures x exposures) and
#' normalise this sum by the total number of mutations in the sample. Then, for each sample, compare its
#' normalised SAD to the normalised SAD of the other samples and check if it is significantly different. In practice,
#' a p-value is computed fitting a gaussian distribution to the other samples.
#' 
#' @param fileout if specified, generate a plot, otherwise no plot is generated
#' @param catalogue original catalogue, channels as rows and samples as columns
#' @param sigs mutational signautures used for fitting, channels as rows, signatures as columns
#' @param exposures exposures/activities of signatures in each sample. Signatures as rows, samples as columns
#' @param pvalue_threshold threshold for statistical significance
#' @keywords unexplained samples
#' @return table of unexplained samples
#' @export
#' @examples
#' res <- SignatureFit_withBootstrap(cat = catalogue,
#'                   signature_data_matrix = cosmic30,
#'                   nboot = 5,
#'                   threshold_percent = 0.1,
#'                   threshold_p.value = 0.1)
#' s_table <- unexplainedSamples(catalogue=catalogue,
#'                   sigs=cosmic30,
#'                   exposures=res$E_median_filtered)
unexplainedSamples <- function(fileout=NULL,catalogue,sigs,exposures,pvalue_threshold=0.01){
  reconstructed <- as.matrix(sigs) %*% as.matrix(exposures)
  
  #look at various metrics, normalised by number of mutation
  norm_SAD <- c()
  for (i in 1:ncol(catalogue)){
    norm_SAD <- c(norm_SAD,sum(abs(reconstructed[,i] - catalogue[,i]))/sum(catalogue[,i]))
  }
  
  pval_norm_SAD <- c()
  for (i in 1:length(norm_SAD)){
    m <- mean(norm_SAD[-i])
    s <- sd(norm_SAD[-i])
    pval_norm_SAD <- c(pval_norm_SAD,1 - pnorm(norm_SAD[i],mean = m,sd = s,lower.tail = TRUE))
  }
  which_significant <- which(pval_norm_SAD<pvalue_threshold)
  
  #plot the outliers
  if(!is.null(fileout)){
    jpeg(filename = fileout,width = 1200,height = 800,res = 200)
    par(mar=c(4,4,5,2))
    plot(norm_SAD,
         col=rgb(0.5,0.5,0.5,0.5),
         pch = 16, ylab = "SAD/n. mutations",xlab = "samples",main = paste0("Normalised Sum of Absolute Deviations (SAD)\n between original and reconstructed catalogue"))
    points(which_significant,norm_SAD[which_significant],col="red",pch = 16)
    legend(x="topleft",legend = c(paste0("significantly higher (p-value<",pvalue_threshold,")")),col = "red",pch = 16,cex = 0.9,bty = "n",inset = c(0,-0.14),xpd = TRUE)
    dev.off()
  }
  #also write the outliers
  res <- data.frame(index=which_significant,sample=colnames(catalogue)[which_significant],normSAD=norm_SAD[which_significant])
  return(res)
}


#' plot_SignatureFit_withBootstrap
#' 
#' Generate plots from the output of SignatureFit_withBootstrap.
#' 
#' @param outdir output directory for the plot
#' @param boostrapFit_res output of SignatureFit_withBootstrap
#' @param type_of_mutations type of mutations: subs, rearr, dnv, or generic
#' @export
plot_SignatureFit_withBootstrap <- function(outdir,
                                          boostrapFit_res,
                                          type_of_mutations){
  
  dir.create(outdir,showWarnings = F,recursive = T)
  
  res <- boostrapFit_res
  
  #function to draw a legend for the heatmap of the correlation matrix
  draw_legend <- function(col,xl,xr,yb,yt){
    par(xpd=TRUE)
    rect(xl,yb,xr,yt)
    rect(
      xl,
      head(seq(yb,yt,(yt-yb)/length(col)),-1),
      xr,
      tail(seq(yb,yt,(yt-yb)/length(col)),-1),
      col=col,border = NA
    )
    text(x = 1.2, y = yt,labels = "1")
    text(x = 1.2, y = (yt-yb)/2,labels = "0")
    text(x = 1.2, y = yb,labels = "-1")
    par(xpd=FALSE)
  }
  
  reconstructed_with_median <- round(as.matrix(res$signature_data_matrix) %*% res$E_median_filtered)
  
  
  #plot and save exposures
  sums_exp <- apply(res$cat, 2, sum)
  exposures <- rbind(res$E_median_filtered,res$unassigned_muts)
  rownames(exposures)[nrow(exposures)] <- "unassigned"
  denominator <- matrix(sums_exp,nrow = nrow(exposures),ncol = ncol(exposures),byrow = TRUE)
  exposuresProp <- (exposures/denominator*100)
  # case of empty catalogues
  exposuresProp[,sums_exp==0] <- 0
  
  file_table_exp <- paste0(outdir,"SigFit_withBootstrap_Exposures_m",res$method,"_bfm",res$bf_method,"_alpha",res$alpha,"_tr",res$threshold_percent,"_p",res$threshold_p.value,".tsv")
  file_plot_exp <- paste0(outdir,"SigFit_withBootstrap_Exposures_m",res$method,"_bfm",res$bf_method,"_alpha",res$alpha,"_tr",res$threshold_percent,"_p",res$threshold_p.value,".jpg")
  file_plot_expProp <- paste0(outdir,"SigFit_withBootstrap_Exposures_m",res$method,"_bfm",res$bf_method,"_alpha",res$alpha,"_tr",res$threshold_percent,"_p",res$threshold_p.value,"_prop.jpg")
  
  plotCosSimMatrix(as.data.frame(exposures),output_file = file_plot_exp,ndigitsafterzero = 0)
  plotCosSimMatrix(as.data.frame(exposuresProp),output_file = file_plot_expProp,ndigitsafterzero = 0)
  write.table(exposures,file = file_table_exp,
              sep = "\t",col.names = TRUE,row.names = TRUE,quote = FALSE)
  
  #provide a series of plots for each sample
  #plot_nrows <- ncol(cat)
  rows_ordered_from_best <- order(res$KLD_samples)
  plot_nrows <- 2
  plot_ncol <- 4
  nplots <- plot_nrows*plot_ncol
  howmanyplots <- ncol(res$cat)
  plotsdir <- paste0(outdir,"SigFit_withBootstrap_Summary_m",res$method,"_bfm",res$bf_method,"_alpha",res$alpha,"_tr",res$threshold_percent,"_p",res$threshold_p.value,"/")
  dir.create(plotsdir,recursive = TRUE,showWarnings = FALSE)
  for(p in 1:howmanyplots){
    current_samples <- p
    jpeg(filename = paste0(plotsdir,"sigfit_bootstrap_",p,"of",howmanyplots,".jpg"),
         width = 640*(plot_ncol),
         height = 480*plot_nrows,
         res = 150)
    par(mfrow=c(plot_nrows,plot_ncol))
    for(i in current_samples){
      fitIsEmpty <- sum(res$E_median_filtered[,i])==0
      unassigned_mut <- sprintf("%.2f",(sum(res$cat[,i,drop=FALSE]) - sum(reconstructed_with_median[,i,drop=FALSE]))/sum(res$cat[,i,drop=FALSE])*100)
      percentdiff <- sprintf("%.2f",sum(abs(res$cat[,i,drop=FALSE] - reconstructed_with_median[,i,drop=FALSE]))/sum(res$cat[,i,drop=FALSE])*100)
      cos_sim <- sprintf("%.2f",cos_sim(res$cat[,i,drop=FALSE],reconstructed_with_median[,i,drop=FALSE]))
      if(type_of_mutations=="subs"){
        #1 original
        plotSubsSignatures(signature_data_matrix = res$cat[,i,drop=FALSE],add_to_titles = "Catalogue",mar=c(6,3,5,2))
        if(!fitIsEmpty){
          #2 reconstructed
          plotSubsSignatures(signature_data_matrix = reconstructed_with_median[,i,drop=FALSE],add_to_titles = "Model",mar=c(6,3,5,2))
          #3 difference
          #plotSubsSignatures(signature_data_matrix = cat[,i,drop=FALSE] - reconstructed_with_median[,i,drop=FALSE],add_to_titles = paste0("Difference, ",percentdiff,"%"),mar=c(6,3,5,2))
          plotSubsSignatures(signature_data_matrix = res$cat[,i,drop=FALSE] - reconstructed_with_median[,i,drop=FALSE],add_to_titles = paste0("Difference\n(CosSim ",cos_sim,", Unassigned ",unassigned_mut,"%)"),mar=c(6,3,5,2),plot_sum = FALSE)
        }
      }else if(type_of_mutations=="rearr"){
        #1 original
        plotRearrSignatures(signature_data_matrix = res$cat[,i,drop=FALSE],add_to_titles = "Catalogue",mar=c(12,3,5,2))
        if(!fitIsEmpty){
          #2 reconstructed
          plotRearrSignatures(signature_data_matrix = reconstructed_with_median[,i,drop=FALSE],add_to_titles = "Model",mar=c(12,3,5,2))
          #3 difference
          #plotRearrSignatures(signature_data_matrix = cat[,i,drop=FALSE] - reconstructed_with_median[,i,drop=FALSE],add_to_titles = paste0("Difference, ",percentdiff,"%"),mar=c(12,3,5,2))
          plotRearrSignatures(signature_data_matrix = res$cat[,i,drop=FALSE] - reconstructed_with_median[,i,drop=FALSE],add_to_titles = paste0("Difference\n(CosSim ",cos_sim,", Unassigned ",unassigned_mut,"%)"),mar=c(12,3,5,2),plot_sum = FALSE)
        }
      }else if(type_of_mutations=="generic"){
        #1 original
        plotGenericSignatures(signature_data_matrix = res$cat[,i,drop=FALSE],add_to_titles = "Catalogue",mar=c(6,3,5,2))
        if(!fitIsEmpty){
          #2 reconstructed
          plotGenericSignatures(signature_data_matrix = reconstructed_with_median[,i,drop=FALSE],add_to_titles = "Model",mar=c(6,3,5,2))
          #3 difference
          #plotGenericSignatures(signature_data_matrix = cat[,i,drop=FALSE] - reconstructed_with_median[,i,drop=FALSE],add_to_titles = paste0("Difference, ",percentdiff,"%"),mar=c(6,3,5,2))
          plotGenericSignatures(signature_data_matrix = res$cat[,i,drop=FALSE] - reconstructed_with_median[,i,drop=FALSE],add_to_titles = paste0("Difference\n(CosSim ",cos_sim,", Unassigned ",unassigned_mut,"%)"),mar=c(6,3,5,2),plot_sum = FALSE)
        }
      }else if(type_of_mutations=="dnv"){
        #1 original
        plotDNVSignatures(signature_data_matrix = res$cat[,i,drop=FALSE],add_to_titles = "Catalogue",mar=c(6,3,5,2))
        if(!fitIsEmpty){
          #2 reconstructed
          plotDNVSignatures(signature_data_matrix = reconstructed_with_median[,i,drop=FALSE],add_to_titles = "Model",mar=c(6,3,5,2))
          #3 difference
          #plotDNVSignatures(signature_data_matrix = cat[,i,drop=FALSE] - reconstructed_with_median[,i,drop=FALSE],add_to_titles = paste0("Difference, ",percentdiff,"%"),mar=c(6,3,5,2))
          plotDNVSignatures(signature_data_matrix = res$cat[,i,drop=FALSE] - reconstructed_with_median[,i,drop=FALSE],add_to_titles = paste0("Difference\n(CosSim ",cos_sim,", Unassigned ",unassigned_mut,"%)"),mar=c(6,3,5,2),plot_sum = FALSE)
        }
      }
      if(!fitIsEmpty){
        #4 bootstraps
        par(mar=c(6,4,5,2))
        boxplot(t(res$samples_list[[i]]),las=3,cex.axes=0.9,
                ylab="n mutations",
                ylim=c(0,max(res$samples_list[[i]])),
                main=paste0("Exposures, of ",colnames(res$E_median_filtered)[i],"\nthreshold=",res$threshold_percent,"%, p-value=",res$threshold_p.value,", n=",res$nboot))
        points(1:length(res$E_median_filtered[,i,drop=FALSE]),res$E_median_filtered[,i,drop=FALSE],col="red")
        abline(h=res$threshold_percent/100*sum(res$cat[,i,drop=FALSE]),col="green")
        legend(x="topleft",legend = c("consensus exposures"),col = "red",pch = 1,cex = 0.9,bty = "n",inset = c(0,-0.14),xpd = TRUE)
        legend(x="topright",legend = c("threshold"),col = "green",lty = 1,cex = 0.9,bty = "n",inset = c(0,-0.14),xpd = TRUE)
        if(ncol(res$signature_data_matrix)>1){
          #5 top correlated signatures
          res.cor <- suppressWarnings(cor(t(res$samples_list[[i]]),method = "spearman"))
          res.cor_triangular <- res.cor
          res.cor_triangular[row(res.cor)+(ncol(res.cor)-col(res.cor))>=ncol(res.cor)] <- 0
          res.cor_triangular_label <- matrix(sprintf("%0.2f",res.cor_triangular),nrow = nrow(res.cor_triangular))
          res.cor_triangular_label[row(res.cor)+(ncol(res.cor)-col(res.cor))>=ncol(res.cor)] <- ""
          par(mar=c(6,8,5,6))
          par(xpd=FALSE)
          col<- colorRampPalette(c("blue", "white", "red"))(51)
          image(res.cor_triangular,col = col,zlim = c(-1,1), axes=F,main="Exposures Correlation (spearman)")
          extrabit <- 1/(ncol(res$signature_data_matrix)-1)/2
          abline(h=seq(0-extrabit,1+extrabit,length.out = ncol(res$signature_data_matrix)+1),col="grey",lty=2)
          abline(v=seq(0-extrabit,1+extrabit,length.out = ncol(res$signature_data_matrix)+1),col="grey",lty=2)
          axis(2,at = seq(0,1,length.out = ncol(res$signature_data_matrix)),labels = colnames(res$signature_data_matrix),las=1,cex.lab=0.8)
          axis(1,at = seq(0,1,length.out = ncol(res$signature_data_matrix)),labels = colnames(res$signature_data_matrix),las=2,cex.lab=0.8)
          draw_legend(col,1.25,1.3,0,1)
          
          #6 some correlation plots
          #pos <- which(max(abs(res.cor_triangular))==abs(res.cor_triangular),arr.ind = TRUE)
          vals <- res.cor_triangular[order(abs(res.cor_triangular),decreasing = TRUE)]
          for (j in 1:(nplots-5)){
            pos <- which(vals[j]==res.cor_triangular,arr.ind = TRUE)
            mainpar <- paste0("Exposures across bootstraps, n=",res$nboot,"\nspearman correlation ",sprintf("%.2f",vals[j]))
            plot(res$samples_list[[i]][pos[1],],res$samples_list[[i]][pos[2],],
                 xlab = colnames(res$signature_data_matrix)[pos[1]],
                 ylab = colnames(res$signature_data_matrix)[pos[2]],
                 # ylim = c(0,max(res$samples_list[[i]][pos[2],])),
                 # xlim = c(0,max(res$samples_list[[i]][pos[1],]))
                 main=mainpar,col="blue",pch = 16)
            
          }
          #sig.pca <- prcomp(t(res$samples_list[[i]]),center = TRUE,scale. = TRUE)
        }
      }
    }
    dev.off()
  }
  # res$E_median_filtered[,i,drop=FALSE]
  
}

