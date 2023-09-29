

#' fitSignaturesWithSignatures
#' 
#' This function attempts to reconstruct a set of signatures (signatures_to_fit) using
#' a second set of signatures (signatures_to use), either just checking the similarity
#' between the signatures of the first and the second set, or using a linear combination
#' of signatures in the second set to model the first set. The typical use of this
#' function would be to check whether a set of extracted signatures can be matched to
#' a set of known signatures, a procedure often called signature assignment or signature
#' identification.
#' 
#' @param signatures_to_fit matrix with the mutational signatures that need to be identified/fitted. Each column has a signature and each row corresponds to a mutational channel. Columns need to sum to 1.
#' @param signatures_to_use matrix of mutational signatures to be used as a reference and that will be fitted/combined to match the signatures_to_fit. Each column has a signature and each row corresponds to a mutational channel. Columns need to sum to 1.
#' @param max_combinations maximum number of signatures_to_use to be used as a linear combination to fit each signatures_to_fit
#' @param sig_min_contrib minimum % of contribution that a signatures_to_use should have to a linear combination when matching a signature in signatures_to_fit
#' @param minsim minimum cosine similarity between a signatures_to_fit and a match, where a match is either a signature in signatures_to_use or a linear combination of signatures in signatures_to_use
#' @param max_from_each_class maximum number of matches to consider in each class, where the class indicates the number of signatures in a linear combination. If there are more matches in the class than the maximum specified, then the matches with the highest cosine similarity will be reported
#' @param nparallel how many parallel cpus to use
#' @return tables with the best linear combinations of signatures in signatures_to_use that match signatures in signatures_to_fit
#' @export
fitSignaturesWithSignatures <- function(signatures_to_fit,
                                        signatures_to_use,
                                        max_combinations = 2,
                                        sig_min_contrib = 0.15,
                                        minsim = 0.75, #min similarity of combination to original
                                        max_from_each_class = 10, #max fits to consider for each class (i.e. number of sigs in a combination)
                                        nparallel = 1){

  doParallel::registerDoParallel(nparallel)
  
  nsigs <- ncol(signatures_to_fit)
  sigs_names <- colnames(signatures_to_fit)
  nsigs_to_use <- ncol(signatures_to_use)

  #get all the combinations
  combinations <- list()
  combinations_class <- c()
  ci <- 1
  for (c in 1:max_combinations){
    combs_table <- combn(1:nsigs_to_use,c)
    for(b in 1:ncol(combs_table)){
      combinations[[ci]] <- combs_table[,b]
      combinations_class <- c(combinations_class,c)
      ci <- ci+1
    }
  }
  
  #go and fit all combinations
  message("[info fitSignaturesWithSignatures] Fitting all the combinations...")
  par_results <- foreach::foreach(i=1:length(combinations)) %dopar% {
    current_sigs <- signatures_to_use[,combinations[[i]],drop=FALSE]
    res <- SignatureFit(signatures_to_fit,current_sigs,doRound = FALSE,verbose = FALSE,showDeprecated = F)
    reconstructed <- as.matrix(current_sigs) %*% res
    cos_sim_fits <- c()
    for(s in 1:nsigs){
      cos_sim_fits <- c(cos_sim_fits,cos_sim(signatures_to_fit[,s],reconstructed[,s]))
    }
    names(cos_sim_fits) <- colnames(signatures_to_fit)
    #return data
    return_res <- list()
    return_res[["fits"]] <- res
    return_res[["cossim"]] <- cos_sim_fits
    return_res
  }
  
  message("[info fitSignaturesWithSignatures] Reorganising the results...")
  #reorganise data
  fit_tables <- list()
  fit_tables_combined <- list()
  combinations_class_list <- list()
  cossim_list <- list()
  for (s in 1:nsigs){
    signatureName <- sigs_names[s]
    fit_tables[[signatureName]] <- NULL
    fit_tables_combined[[signatureName]] <- NULL
    combinations_class_list[[signatureName]] <- c()
    cossim_list[[signatureName]] <- c()
    for (i in 1:length(combinations)){
      #discard the clearly poor fits here, also discard if one of the sigs doesn't contribute much
      if(par_results[[i]]$cossim[signatureName] >= minsim & all(par_results[[i]]$fits[,signatureName]>=sig_min_contrib)){ 
        newrow <- vector(mode = "numeric",length = nsigs_to_use)
        names(newrow) <- colnames(signatures_to_use)
        newrow[combinations[[i]]] <- par_results[[i]]$fits[,signatureName]
        fit_tables[[signatureName]] <- rbind(fit_tables[[signatureName]],newrow)
        combinations_class_list[[signatureName]] <- c(combinations_class_list[[signatureName]],combinations_class[i])
        cossim_list[[signatureName]] <- c(cossim_list[[signatureName]],par_results[[i]]$cossim[signatureName])
      }
    }

    #get the best of each class
    best_table <- NULL
    best_class <- c()
    best_cossim <- c()
    for (c in 1:max_combinations){
      tmp_slice <- fit_tables[[signatureName]][combinations_class_list[[signatureName]]==c,,drop=FALSE]
      if(!is.null(tmp_slice)){
        if (nrow(tmp_slice)>0){
          tmp_cossim <- cossim_list[[signatureName]][combinations_class_list[[signatureName]]==c]
          row_selection <- order(tmp_cossim,decreasing = TRUE)[1:min(nrow(tmp_slice),max_from_each_class)]
          best_cossim <- c(best_cossim,tmp_cossim[row_selection])
          best_class <- c(best_class,rep(c,length(row_selection)))
          best_table <- rbind(best_table,tmp_slice[row_selection,,drop=FALSE])
        }
      }
    }
    #replace
    fit_tables[[signatureName]] <- best_table
    combinations_class_list[[signatureName]] <- best_class
    cossim_list[[signatureName]] <- best_cossim
    
    if(nrow(fit_tables[[signatureName]])>0){
      fit_tables_combined[[signatureName]] <- cbind(data.frame(nsigs=combinations_class_list[[signatureName]],
                                                               cossim=cossim_list[[signatureName]],
                                                               stringsAsFactors = F),
                                                    fit_tables[[signatureName]])
      rownames(fit_tables_combined[[signatureName]]) <- 1:nrow(fit_tables_combined[[signatureName]])
    }
  }
  
  final_res <- list()
  final_res[["fit_tables"]] <- fit_tables
  final_res[["combinations_class_list"]] <- combinations_class_list
  final_res[["cossim_list"]] <- cossim_list
  final_res[["fit_tables_combined"]] <- fit_tables_combined
  return(final_res)
  
}

