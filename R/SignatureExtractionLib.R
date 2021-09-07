#Fork of SigHunter from Sandro Morganella 2017
#Andrea Degasperi, ad923@cam.ac.uk



#' Mutational Signatures Extraction
#' 
#' Perform signature extraction, by applying NMF to the input matrix. Multiple NMF runs and bootstrapping is used for robustness, followed by clustering of the solutions. A range of number of signatures to be used is required.
#' 
#' @param cat matrix with samples as columns and channels as rows
#' @param matrix_of_fixed_signatures matrix with known signatures as columns and channels as rows. Used for partial extraction with NNLM package, with Lee KLD (brunet) only. If NULL, NMF package is used instead and different nmf methods can be used.
#' @param outFilePath path were the extraction output files should go. Remember to add "/" at the end of the path
#' @param blacklist list of samples (column names) to ignore
#' @param nrepeats how many runs for each bootstrap (if filterBestOfEachBootstrap=TRUE with default params, only at most 10 runs within 0.1 percent of best will be considered, so nrepeats should be at least 10)
#' @param nboots how many bootstrapped catalogues to use
#' @param clusteringMethod choose among {"HC","PAM","MC"}, hierarchical clustering (HC), partitioning around the medoids (PAM) and  matched clustering (MC)
#' @param completeLinkageFlag if clusteringMethod="HC", use complete linkage instead of default average linkage
#' @param useMaxMatching if clusteringMethod="MC", use the assignment problem algorithm (match with max similarity) instead of the stable matching algorithm (any stable match)
#' @param filterBestOfEachBootstrap if TRUE only at most filterBest_nmaxtokeep of the nrepeats runs that are within filterBest_RTOL*best from the best are kept
#' @param filterBest_RTOL realtive tolerace from best fit to consider a run as good as the best, RTOL=0.001 is recommended
#' @param filterBest_nmaxtokeep max number of runs that should be kept that are within the relative tolerance from the best
#' @param nparallel how many processing units to use
#' @param nsig list of number of signatures to try
#' @param mut_thr threshold of mutations to remove empty/almost empty rows and columns
#' @param type_of_extraction choose among {"subs","rearr","generic","dnv"}
#' @param project give a name to your project
#' @param parallel set to TRUE to use parallel computation (Recommended)
#' @param nmfmethod choose among {"brunet","lee","nsNMF"}, this choice will be passed to the NMF::nmf function
#' @param removeDuplicatesInCatalogue remove 0.99 cos sim similar samples
#' @param normaliseCatalogue scale samples to sum to 1
#' @param plotCatalogue also plot the catalogue, this may crash the library if the catalogue is too big, should work up to ~300 samples
#' @param plotResultsFromAllClusteringMethods if TRUE, all clustering methods are used and results are reported and plotted for all of them. If FALSE, only the requested clustering is reported
#' @return result files will be available in the outFilePath directory
#' @keywords signature extraction
#' @export
#' @examples
#'   n_row <- 96
#'   n_col <- 50
#'   rnd_matrix <- round(matrix(runif(n_row*n_col,min = 0,max = 50),nrow = n_row,ncol = n_col))
#'   colnames(rnd_matrix) <- paste0("C",1:n_col)
#'   row.names(rnd_matrix) <- paste0("R",1:n_row)
#'   SignatureExtraction(cat = rnd_matrix,
#'                       outFilePath = paste0("extraction_test_subs/"),
#'                       nrepeats = 10,
#'                       nboots = 2,
#'                       nparallel = 2,
#'                       nsig = 2:3,
#'                       mut_thr = 0,
#'                       type_of_extraction = "subs",
#'                       project = "test",
#'                       parallel = TRUE,
#'                       nmfmethod = "brunet")
SignatureExtraction <- function(cat, #matrix with samples as columns and channels as rows.
                                outFilePath, #path were the extraction output files should go. Remember to add "/" at the end of the path
                                matrix_of_fixed_signatures=NULL, #matrix with known signatures as columns and channels as rows.
                                blacklist=c(), #list of samples (column names) to ignore
                                nrepeats=10, #how many runs for each bootstrap (if filterBestOfEachBootstrap=TRUE with default params, only at most 10 runs within 0.1% of best will be considered, so nrepeats should be at least 10)
                                nboots=20, #how many bootstrapped catalogues to use
                                clusteringMethod="MC", #choose among {"HC","PAM","MC"}, hierarchical clustering (HC), partitioning around the medoids (PAM) and  matched clustering (MC)
                                completeLinkageFlag=FALSE, #if clusteringMethod="HC", use complete linkage instead of default average linkage
                                useMaxMatching=TRUE, #if clusteringMethod="MC", use the assignment problem algorithm (match with max similarity) instead of the stable matching algorithm (any stable match)
                                filterBestOfEachBootstrap=TRUE, #if true only at most filterBest_nmaxtokeep of the nrepeats runs that are within filterBest_RTOL*best from the best are kept
                                filterBest_RTOL=0.001, #realtive tolerace from best fit to consider a run as good as the best
                                filterBest_nmaxtokeep=10, #max number of runs that should be kept that are within the relative tolerance from the best
                                nparallel=1, # how many processing units to use
                                nsig=c(3:15), # range of number of signatures to try
                                mut_thr=0, # threshold of mutations to remove empty/almost empty rows and columns
                                type_of_extraction="subs", #choose among {"subs","rearr","generic","dnv"}
                                project="extraction", # give a name to your project
                                parallel=FALSE, # set to TRUE to use parallel computation (Recommended)
                                nmfmethod="brunet", #choose among {"brunet","lee","nsNMF"}
                                removeDuplicatesInCatalogue = FALSE, #remove 0.99 cos sim similar samples
                                removeDuplicatesThreshold = 0.98,
                                normaliseCatalogue = FALSE, # scale samples to sum to 1
                                plotCatalogue = FALSE, #also plot the catalogue, this may crash the library if the catalogue is too big, should work up to ~300 samples
                                plotResultsFromAllClusteringMethods=TRUE){ #if TRUE, all clustering methods are used and results are reported and plotted for all of them. If FALSE, only the requested clustering is reported
  
  library(NMF)
  
  #tmp for debug
  if (completeLinkageFlag) tmp_outFilePath <- outFilePath
  
  #remove blacklisted samples
  if(length(blacklist)>0){
    cat <- cat[,setdiff(colnames(cat),blacklist)]
  }
  
  group <- project
  
  doParallel::registerDoParallel(nparallel)
  dir.create(outFilePath,showWarnings = FALSE,recursive = TRUE)
  
  message("\n------------- START COMPUTATION ------------\n")
  
  ## Read the catalogue
  message("[STEP 1]: Reading and Preprocessing the Catalogue...", appendLF=F)
  if (type_of_extraction=="subs") cat <- sortCatalogue(cat)
  if(removeDuplicatesInCatalogue){
    cat <- removeSimilarCatalogueSamples(cat,cosSimThreshold = removeDuplicatesThreshold)
  }
  #Save catalogue used
  nsamples <- ncol(cat)
  if(plotCatalogue){
    cat_file <- paste0(outFilePath,"CatalogueUsedAfterPreprocessing_plot_",project,".pdf")
    if (type_of_extraction == "subs"){
      plotSubsSignatures(cat,cat_file,plot_sum = TRUE,overall_title = paste0("Catalogue - ",project," (",nsamples," samples)"))
    }else if (type_of_extraction == "rearr"){
      plotRearrSignatures(cat,cat_file,plot_sum = TRUE,overall_title = paste0("Catalogue - ",project," (",nsamples," samples)"))
    }else if (type_of_extraction == "generic"){
      plotGenericSignatures(cat,cat_file,plot_sum = TRUE,overall_title = paste0("Catalogue - ",project," (",nsamples," samples)"))
    }else if (type_of_extraction == "dnv"){
      plotDNVSignatures(cat,cat_file,plot_sum = TRUE,overall_title = paste0("Catalogue - ",project," (",nsamples," samples)"))
    }
  }
  cat_file <- paste0(outFilePath,"CatalogueUsedAfterPreprocessing_plot_",group,".txt")
  write.table(cat,file = cat_file,
              sep = "\t",quote = FALSE,row.names = TRUE,col.names = TRUE)
  
  nrow_cat <- nrow(cat)
  
  #remove channels if necessary
  all_rows_cat <- cat
  #cat <-  preprocessCatalgue(cat, mut_thr)
  channelsRemoved <- FALSE
  if(nrow(all_rows_cat)>nrow(cat)) channelsRemoved <- TRUE
  
  #normalised catalogue for computing replicates later
  ncat <- cat
  if(normaliseCatalogue) ncat <- normaliseSamples(cat)
  
  if(sum(nsig>ncol(cat))>0){
    nsig <- nsig[nsig[length(nsig)]<=ncol(cat)]
  }
  
  message("DONE")
  message("\t> ", ncol(cat), " Samples and ",nrow_cat, " Mutations Correclty Processed")
  

  
  if(mut_thr>0){
    message("\t> ", nrow_cat -nrow(cat), " Mutation(s) Channel(s) Removed (<", mut_thr, ")")
  }
  
  nmuts <- apply(cat, 2, sum)
  
  sample_names <- colnames(cat)
  
  ###### Find Signature ######
  message("\n[STEP 2]: Finding Signatures")
  
  ## Iterate Until Convergence
  complete <- FALSE
  round <- 1
  solutions <- matrix(0, 0, 4)
  colnames(solutions) <- c("Round", "Optim_Nsig", "Estimate", "NumUncorrelated") 
  
  prev_num_uncorr <- ncol(cat)
  
  while(complete==FALSE){
    message("\t> Running round ", round)
    
    ## Create Directory
    outDir <- paste(outFilePath, "round_", round, "/", sep="")
    cmd_res <- system(paste("mkdir -p", outDir))
    
    
    err <- c()
    cos_sim <- c()
    average_corr_smpls <- c()
    sil <- c()
    
    #additional metrics
    ave.RMSE <- c()
    sd.RMSE <- c()
    ave.KLD <- c()
    sd.KLD <- c()
    ave.RMSE.orig <- c()
    sd.RMSE.orig <- c()
    ave.KLD.orig <- c()
    sd.KLD.orig <- c()
    # ave.CosSim.hclust <- c()
    # ave.CosSim.PAM <- c()
    ave.SilWid.hclust <- c()
    ave.SilWid.PAM <- c()
    ave.SilWid.MC <- c()
    cophenetic.corr.hclust <- c()
    # proportion.tooSimilar.Signatures <- c()
    min.MinWCCS.hclust <- c()
    min.MinWCCS.PAM <- c()
    min.MinWCCS.MC <- c()
    max.MaxBCCS.hclust <- c()
    max.MaxBCCS.PAM <- c()
    max.MaxBCCS.MC <- c()
    mmcs_pam <- c()
    mmcs_hclust <- c()
    mmcs_MC <- c()
    
    ## Run Hunter on the specified range of signatures
    for(ns in nsig){
      
      #tmp for debug
      if (completeLinkageFlag) outFilePath <- tmp_outFilePath
      
      outNsDir <- paste(outDir, "sig_", ns, "/", sep="")
      cmd_res <- system(paste("mkdir -p ", outNsDir))
      
      message("\n> Running for ", ns,  " Signatures (", nboots, " bootstraps)")
      strt<-Sys.time()
      
      #---------Bootstraps start here
      
      #Define smoothing matrix S in case nsNMF is used
      if (nmfmethod=="nsNMF") S <- diag(ns)*0.5 + matrix(1,nrow=ns,ncol=ns)*0.5/ns
      
      #bootstraps file:
      bootstraps_file <- paste0(outFilePath,"bootstraps_",group,"_ns",ns,"_nboots",nboots,".Rdata")
      if (file.exists(bootstraps_file)){
        load(bootstraps_file)
        message("bootstraps file loaded")
      }else{
        ## Collect the solutions
        e_boot <- matrix(0, ncol(cat), 0)
        p_boot <- matrix(0, nrow(cat), 0)
        err_boot <- c()
        boot_tracker <- c()
        boot_tracker_e <- c()
        #boots_list <- list()
        boot_cat <- list()
        ## Run NMF on nboots bootsrapped catalogues
        if(parallel==TRUE){ ## Parallel
          #boots_list_tmp <- list()
          nseq <- ceiling(nboots/nparallel)
          #first get and store catalogues
          for (i in 1:(nseq*nparallel)){
            #boot_cat[[i]] <- apply(cat, 2, generateRandMuts)
            boot_cat[[i]] <- generateRandMuts(cat)
            if(normaliseCatalogue) boot_cat[[i]] <- normaliseSamples(boot_cat[[i]])
          }
          for(tt in 1:nseq){
            message(".",  appendLF=F)	
            boots_list <- foreach::foreach(i=1:nparallel) %dopar%{						
              rnd_cat <- boot_cat[[(tt-1)*nparallel + i]]			
              rnd_cat <- preprocessCatalgue(rnd_cat, mut_thr) #remove channels with 0 mutations
              if(!is.null(matrix_of_fixed_signatures)){
                nnmf.res <- list()
                for (i in 1:nrepeats){
                  nnmf.res[[i]] <- NNLM::nnmf(as.matrix(rnd_cat),
                                        init = list(W0 = as.matrix(matrix_of_fixed_signatures[rownames(rnd_cat),,drop=F])),
                                        loss = "mkl",
                                        method = "lee",
                                        k = ns,
                                        max.iter = 10000,check.k = FALSE)
                }
                nnmf.res
              }else{
                NMF::nmf(rnd_cat, rank=ns,nrun=nrepeats,method = nmfmethod,.options="k-p")
              }
            }
            #Extract data already to save memory
            for(i in 1:length(boots_list)){
              nmf_res <- boots_list[[i]]
              # get residuals and best residuals
              if(!is.null(matrix_of_fixed_signatures)){
                residuals_list <- c()
                for (j in 1:length(nmf_res)){
                  #avoid numerical error (minimum cannot be less than 0)
                  residuals_list <- c(residuals_list,max(0,nmf_res[[j]]$mkl[length(nmf_res[[j]]$mkl)]))
                }
                best_residual <- min(residuals_list)
              }else{
                best_residual <- NMF::residuals(nmf_res)
                residuals_list <- c()
                if(length(nmf_res)==1){
                  residuals_list <- NMF::residuals(nmf_res)
                }else{
                  for (j in 1:length(nmf_res@.Data)){
                    #avoid numerical error (minimum cannot be less than 0)
                    residuals_list <- c(residuals_list,max(0,NMF::residuals(nmf_res@.Data[[j]])))
                  }
                }
              }
              # decide what to keep
              if(filterBestOfEachBootstrap){
                #filter the best runs
                runsToChooseFrom <- which(best_residual*(1+filterBest_RTOL)>=residuals_list)
                #take at most filterBest_nmaxtokeep
                if (length(runsToChooseFrom)>filterBest_nmaxtokeep){
                  runsToChooseFrom <- sample(runsToChooseFrom,filterBest_nmaxtokeep)
                }
              }else{
                #just take all the runs
                runsToChooseFrom <- 1:length(nmf_res)
              }
              
              countRuns <- 1
              for (j in 1:length(nmf_res)){
                #keep at most 10 repeats that have residuals close to the best (within 0.1% more than residual)
                #if (best_residual*1.001>residuals_list[j] & countRuns <= 10) {
                if (j %in% runsToChooseFrom) {
                  if(length(nmf_res)==1){
                    #here is where I fix channels by adding back the channels I removed
                    if(!is.null(matrix_of_fixed_signatures)){
                      bbb <- nmf_res[[j]]$W[,1:ns,drop=F]
                    }else{
                      bbb <- NMF::basis(nmf_res)
                    }
                    coln <- paste0("b",i,"r",j,"s",1:ns)
                    p_boot_tmp <- matrix(0,nrow = nrow(cat),ncol = ns,dimnames = list(rownames(cat),coln))
                    p_boot_tmp[rownames(bbb),] <- bbb
                    if (nmfmethod=="nsNMF"){
                      p_boot <- cbind(p_boot , p_boot_tmp %*% S)
                    }else{
                      p_boot <- cbind(p_boot , p_boot_tmp)
                    }
                    if(!is.null(matrix_of_fixed_signatures)){
                      e_boot <- cbind(e_boot , t(nmf_res[[j]]$H))
                      err_boot <- c(err_boot,nmf_res[[j]]$mkl[length(nmf_res[[j]]$mkl)])
                    }else{
                      e_boot <- cbind(e_boot , t(NMF::coef(nmf_res)))
                      err_boot <- c(err_boot,NMF::residuals(nmf_res))
                    }
                  }else{
                    #here is where I fix channels by adding back the channels I removed
                    if(!is.null(matrix_of_fixed_signatures)){
                      bbb <- nmf_res[[j]]$W[,1:ns,drop=F]
                    }else{
                      bbb <- NMF::basis(nmf_res@.Data[[j]])
                    }
                    coln <- paste0("b",i,"r",j,"s",1:ns)
                    p_boot_tmp <- matrix(0,nrow = nrow(cat),ncol = ns,dimnames = list(rownames(cat),coln))
                    p_boot_tmp[rownames(bbb),] <- bbb
                    if (nmfmethod=="nsNMF"){
                      p_boot <- cbind(p_boot , p_boot_tmp %*% S)
                    }else{
                      p_boot <- cbind(p_boot , p_boot_tmp)
                    }
                    if(!is.null(matrix_of_fixed_signatures)){
                      e_boot <- cbind(e_boot , t(nmf_res[[j]]$H))
                      err_boot <- c(err_boot,nmf_res[[j]]$mkl[length(nmf_res[[j]]$mkl)])
                    }else{
                      e_boot <- cbind(e_boot , t(NMF::coef(nmf_res@.Data[[j]])))
                      err_boot <- c(err_boot,NMF::residuals(nmf_res@.Data[[j]]))
                    }
                    
                  }

                  boot_tracker <- c(boot_tracker,rep((tt-1)*nparallel + i,ns))
                  if(!is.null(matrix_of_fixed_signatures)){
                    boot_tracker_e <- c(boot_tracker_e,rep((tt-1)*nparallel + i,ns+ncol(matrix_of_fixed_signatures)))
                  }else{
                    boot_tracker_e <- c(boot_tracker_e,rep((tt-1)*nparallel + i,ns))
                  }
                  countRuns <- countRuns + 1
                }
              }
              rm(nmf_res)
              gc()
            }
            rm(boots_list)
            gc()
          }
          #Now you need to cut so that you have only nboots and not more
          selection_of_boots <- boot_tracker <= nboots
          selection_of_boots_e <- boot_tracker_e <= nboots
          selection_of_boots_runs <- boot_tracker[seq(1,length(boot_tracker),ns)] <= nboots
          p_boot <- p_boot[,selection_of_boots]
          e_boot <- e_boot[,selection_of_boots_e]
          err_boot <- err_boot[selection_of_boots_runs]
          boot_tracker <- boot_tracker[selection_of_boots]
          # boot_pos <- 1
          # for (s in 1:nseq){
          #   for (p in 1:length(boots_list_tmp[[s]])){
          #     boots_list[[boot_pos]] <- boots_list_tmp[[s]][[p]]
          #     boot_pos <- boot_pos + 1
          #   }
          # }
          #boots_list <- unlist(boots_list_tmp)
        }else{ ## No Parallel	
          for(i in 1:nboots){
            message(".",  appendLF=F)			
            boot_cat[[i]] <- generateRandMuts(cat)
            if(normaliseCatalogue) boot_cat[[i]] <- normaliseSamples(boot_cat[[i]])
            rnd_cat <- preprocessCatalgue(boot_cat[[i]], mut_thr) #remove channels with 0 mutations
            if(!is.null(matrix_of_fixed_signatures)){
              nmf_res <- list()
              for (j in 1:nrepeats){
                nmf_res[[j]] <- NNLM::nnmf(as.matrix(rnd_cat),
                                      init = list(W0 = as.matrix(matrix_of_fixed_signatures[row.names(rnd_cat),,drop=F])),
                                      loss = "mkl",
                                      method = "lee",
                                      k = ns,
                                      max.iter = 10000,check.k = FALSE)
              }
            }else{
              nmf_res <- NMF::nmf(rnd_cat, rank=ns,nrun=nrepeats,method = nmfmethod,.options="k-p")
            }
            #extract already to save memory
            if(!is.null(matrix_of_fixed_signatures)){
              residuals_list <- c()
              for (j in 1:length(nmf_res)){
                #avoid numerical error (minimum cannot be less than 0)
                residuals_list <- c(residuals_list,max(0,nmf_res[[j]]$mkl[length(nmf_res[[j]]$mkl)]))
              }
              best_residual <- min(residuals_list)
            }else{
              best_residual <- NMF::residuals(nmf_res)
              residuals_list <- c()
              if(length(nmf_res)==1){
                residuals_list <- NMF::residuals(nmf_res)
              }else{
                for (j in 1:length(nmf_res@.Data)){
                  #avoid numerical error (minimum cannot be less than 0)
                  residuals_list <- c(residuals_list,max(0,NMF::residuals(nmf_res@.Data[[j]])))
                }
              }
            }
            
            if(filterBestOfEachBootstrap){
              #filter the best runs
              runsToChooseFrom <- which(best_residual*(1+filterBest_RTOL)>=residuals_list)
              #take at most filterBest_nmaxtokeep
              if (length(runsToChooseFrom)>filterBest_nmaxtokeep){
                runsToChooseFrom <- sample(runsToChooseFrom,filterBest_nmaxtokeep)
              }
            }else{
              #just take all the runs
              runsToChooseFrom <- 1:length(nmf_res)
            }
            countRuns <- 1
            for (j in 1:length(nmf_res)){
              #keep at most 10 repeats that have residuals close to the best (within 0.1% more than residual)
              #if (best_residual*1.001>residuals_list[j] & countRuns <= 10) {
              if (j %in% runsToChooseFrom) {
                #here is where I fix channels by adding back the channels I removed
                if(length(nmf_res)==1){
                  if(!is.null(matrix_of_fixed_signatures)){
                    bbb <- nmf_res[[j]]$W[,1:ns,drop=F]
                  }else{
                    bbb <- NMF::basis(nmf_res)
                  }
                  coln <- paste0("b",i,"r",j,"s",1:ns)
                  p_boot_tmp <- matrix(0,nrow = nrow(cat),ncol = ns,dimnames = list(rownames(cat),coln))
                  p_boot_tmp[rownames(bbb),] <- bbb
                  if (nmfmethod=="nsNMF"){
                    p_boot <- cbind(p_boot , p_boot_tmp %*% S)
                  }else{
                    p_boot <- cbind(p_boot , p_boot_tmp)
                  }
                  if(!is.null(matrix_of_fixed_signatures)){
                    e_boot <- cbind(e_boot , t(nmf_res[[j]]$H))
                    err_boot <- c(err_boot,nmf_res[[j]]$mkl[length(nmf_res[[j]]$mkl)])
                  }else{
                    e_boot <- cbind(e_boot , t(NMF::coef(nmf_res)))
                    err_boot <- c(err_boot,NMF::residuals(nmf_res))
                  }
                }else{
                  #here is where I fix channels by adding back the channels I removed
                  if(!is.null(matrix_of_fixed_signatures)){
                    bbb <- nmf_res[[j]]$W[,1:ns,drop=F]
                  }else{
                    bbb <- NMF::basis(nmf_res@.Data[[j]])
                  }
                  coln <- paste0("b",i,"r",j,"s",1:ns)
                  p_boot_tmp <- matrix(0,nrow = nrow(cat),ncol = ns,dimnames = list(rownames(cat),coln))
                  p_boot_tmp[rownames(bbb),] <- bbb
                  if (nmfmethod=="nsNMF"){
                    p_boot <- cbind(p_boot , p_boot_tmp %*% S)
                  }else{
                    p_boot <- cbind(p_boot , p_boot_tmp)
                  }
                  if(!is.null(matrix_of_fixed_signatures)){
                    e_boot <- cbind(e_boot , t(nmf_res[[j]]$H))
                    err_boot <- c(err_boot,nmf_res[[j]]$mkl[length(nmf_res[[j]]$mkl)])
                  }else{
                    e_boot <- cbind(e_boot , t(NMF::coef(nmf_res@.Data[[j]])))
                    err_boot <- c(err_boot,NMF::residuals(nmf_res@.Data[[j]]))
                  }
                  
                }
                boot_tracker <- c(boot_tracker,rep(i,ns))
                boot_tracker_e <- c(boot_tracker_e,rep(i,ns+ncol(matrix_of_fixed_signatures)))
                countRuns <- countRuns + 1
              }
            }
            rm(nmf_res)
            gc()
          }	
        }
        #how many solutions were saved?
        saved_nmf_runs <- ncol(p_boot)/ns
        save(file = bootstraps_file,ns,nboots,nrepeats,e_boot,p_boot,err_boot,cat,all_rows_cat,saved_nmf_runs,boot_tracker,boot_tracker_e,boot_cat)
      }
      
      #now add back the missing channel and reset cat to all channels
      
      # if(channelsRemoved){
      #   lmissing <- setdiff(rownames(all_rows_cat),rownames(cat))
      #   nmissing <- length(lmissing)
      #   newrows <- matrix(0,nrow = nmissing,ncol = ncol(p_boot))
      #   colnames(newrows) <- colnames(p_boot)
      #   rownames(newrows) <- lmissing
      #   p_boot <- rbind(p_boot,newrows)[rownames(all_rows_cat),]
      #   #cat <- all_rows_cat
      # }
      
      # ## Compute the average silhouette grouping all the computed solutions in ns clusters
      #sil <- c(sil, summary(silhouette(pam(p_boot, ns)))$avg.width)
      #above dissimilarity is not based on cosine similarity
      colnames_p <- c()
      sigs_p <- 1:ns
      for (i in 1:saved_nmf_runs){
        colnames_p <- c(colnames_p,paste(rep("sig",ns),sigs_p,rep("_run",ns),rep(i,ns),sep = ""))
      }
      colnames(p_boot) <- colnames_p
      names(boot_tracker) <- colnames_p
      #Load distance matrix if it already exists
      distMatrix_file <- paste0(outFilePath,"distMatrix_",group,"_ns",ns,"_nboots",nboots,".Rdata")
      if (file.exists(distMatrix_file)){
        load(distMatrix_file)
        message("distance matrix file loaded")
      }else{
        distMatrix <- 1 - computeCorrelation_parallel(p_boot,nparallel = nparallel,parallel = TRUE)
        save(file = distMatrix_file,distMatrix)
      }
      
      #change prefix (debug)
      if (completeLinkageFlag) outNsDir <- paste0(outNsDir,"completeLinkage")
      
      if(ns>1 & saved_nmf_runs>1){
        #Hierarchical clustering
        if (completeLinkageFlag) {
          fit_clust <- stats::hclust(as.dist(distMatrix), method="complete") 
        }else{
          fit_clust <- stats::hclust(as.dist(distMatrix), method="average") 
        }
        #Hierarchical clustering partitioning
        cut_res <- stats::cutree(fit_clust,k = ns)
        #PAM partitioning
        clustering_result <- cluster::pam(as.dist(distMatrix), ns)
        cut_res.PAM <- clustering_result$clustering
        #Matched Clustering
        mc_file <- paste0(outNsDir,"matchedClustering_",group,"_ns",ns,"_nboots",nboots,".Rdata")
        if (file.exists(mc_file)){
          load(mc_file)
          message("matched clustering result loaded from file")
        }else{
          cut_res_MC <- matchedClustering(distMatrix,ns,maxMatch = useMaxMatching,parallel=TRUE,nparallel=nparallel)
          save(file = mc_file,cut_res_MC)
        }
        #get medoids
        medoids_hclust <- findMedoidsHclust(distMatrix,cut_res)
        medoids_MC <- findMedoidsHclust(distMatrix,cut_res_MC)
        medoids_PAM <- clustering_result$medoids
        
      }else if(ns==1 & saved_nmf_runs>1){
        cut_res <- rep(1,ncol(distMatrix))
        names(cut_res) <- colnames(distMatrix)
        cut_res.PAM <- cut_res
        cut_res_MC <- cut_res
        medoids_hclust <- findMedoidsHclust(distMatrix,cut_res)
        medoids_MC <- medoids_hclust
        medoids_PAM <- medoids_hclust
      }else if(ns>1 & saved_nmf_runs==1){
        cut_res <- seq(1,ncol(distMatrix))
        names(cut_res) <- colnames(distMatrix)
        cut_res.PAM <- cut_res
        cut_res_MC <- cut_res
        medoids_hclust <- findMedoidsHclust(distMatrix,cut_res)
        medoids_MC <- medoids_hclust
        medoids_PAM <- medoids_hclust
      }else if(ns==1 & saved_nmf_runs==1){
        cut_res <- 1
        names(cut_res) <- colnames(distMatrix)
        cut_res.PAM <- cut_res
        cut_res_MC <- cut_res
        medoids_hclust <- colnames(distMatrix)
        medoids_MC <- medoids_hclust
        medoids_PAM <- medoids_hclust
      }
      
      norm_p_boot <- p_boot/matrix(data=rep(apply(p_boot,2,sum),nrow(p_boot)),nrow = nrow(p_boot),byrow = TRUE)
      
      #mean and sd
      #set up mean and sd of signatures
      mean_signatures <- list()
      sd_signatures <- list()
      for (cm in c("HC","PAM","MC")){
        mean_signatures[[cm]] <- matrix(NA,nrow = nrow(p_boot),ncol = ns)
        sd_signatures[[cm]] <- matrix(NA,nrow = nrow(p_boot),ncol = ns)
        colnames(mean_signatures[[cm]]) <- paste0("S",1:ns)
        colnames(sd_signatures[[cm]]) <- paste0("S",1:ns)
        row.names(mean_signatures[[cm]]) <- row.names(p_boot)
        row.names(sd_signatures[[cm]]) <- row.names(p_boot)
      }
      
      for (cm in c("HC","PAM","MC")){
        if(cm=="HC"){
          partitions <- cut_res
        }else if(cm=="PAM"){
          partitions <- cut_res.PAM
        }else if(cm=="MC"){
          partitions <- cut_res_MC
        }
        for(nsi in 1:ns){
          pboot_dim <- dim(norm_p_boot[,partitions==nsi,drop=FALSE])
          if(pboot_dim[2]==1){
            mean_signatures[[cm]][,nsi] <- norm_p_boot[,partitions==nsi]
            sd_signatures[[cm]][,nsi] <- 0
          }else{
            mean_signatures[[cm]][,nsi] <- apply(norm_p_boot[,partitions==nsi],1,mean)
            sd_signatures[[cm]][,nsi] <- apply(norm_p_boot[,partitions==nsi],1,sd)
            sd_signatures[[cm]][is.na(sd_signatures[[cm]][,nsi]),nsi] <- 0
          }
        }
        if(!is.null(matrix_of_fixed_signatures)){
          mean_signatures[[cm]] <- cbind(matrix_of_fixed_signatures,mean_signatures[[cm]])
          tmpmatrix <- matrix_of_fixed_signatures
          tmpmatrix[,] <- 0
          sd_signatures[[cm]] <- cbind(tmpmatrix,sd_signatures[[cm]])
        }
      }
      
      if(ns>1 & saved_nmf_runs>1){
        #Compute Silhouettes
        sil_hclust <- summary(cluster::silhouette(cut_res,as.dist(distMatrix)))
        sil_pam <- summary(cluster::silhouette(cut_res,as.dist(distMatrix)))
        sil_MC <- summary(cluster::silhouette(cut_res_MC,as.dist(distMatrix)))
        plotWithinClusterSilWidth(sil_hclust,sil_pam,sil_MC,outNsDir,group,ns,nboots)
        
        #compute within cluster cos similiraty distance
        cosSimPAM <- withinClusterCosSim(clustering_result$clustering,distMatrix,parallel = TRUE)
        cosSimHClust <- withinClusterCosSim(cut_res,distMatrix,parallel = TRUE)
        cosSimMC <- withinClusterCosSim(cut_res_MC,distMatrix,parallel = TRUE)
        plotWithinClusterCosSim(cosSimHClust,cosSimPAM,cosSimMC,outNsDir,group,ns,nboots)
        
        #compute cophenetic correlation
        coph <- cor(stats::cophenetic(fit_clust),as.dist(distMatrix))
        
        #compute the proportion of the solutions that contain signatures that are too similar
        #propTooSimilar <- computePropTooSimilar(distMatrix,saved_nmf_runs,ns)
        
        #Save metrics
        #ave.CosSim.hclust <- c(ave.CosSim.hclust,mean(cosSimHClust))
        #ave.CosSim.PAM <- c(ave.CosSim.PAM,mean(cosSimPAM))
        ave.SilWid.hclust <- c(ave.SilWid.hclust,sil_hclust$avg.width)
        ave.SilWid.PAM <- c(ave.SilWid.PAM,sil_pam$avg.width)
        ave.SilWid.MC <- c(ave.SilWid.MC,sil_MC$avg.width)
        cophenetic.corr.hclust <- c(cophenetic.corr.hclust,coph)
        #proportion.tooSimilar.Signatures <- c(proportion.tooSimilar.Signatures,propTooSimilar)
        
        #save metrics
        additional_perf_file <- paste0(outNsDir,"Sigs_WithinClusterPerf_",group,"_ns",ns,"_nboots",nboots,".rData")
        save(file = additional_perf_file,sil_pam,sil_hclust,sil_MC,cosSimPAM,cosSimHClust,cosSimMC)
        #plotHierarchicalCluster(fit_clust,outNsDir,group,ns,nboots)
        
        #Max Medoids Cosine Similarity
        mmcs_pam <- c(mmcs_pam,max(medoids_cosSimMatrix(p_boot,medoids_PAM) - diag(ns)))
        mmcs_hclust <- c(mmcs_hclust,max(medoids_cosSimMatrix(p_boot,medoids_hclust) - diag(ns)))
        mmcs_MC <- c(mmcs_MC,max(medoids_cosSimMatrix(p_boot,medoids_MC) - diag(ns)))
        
      }else{
        ave.SilWid.hclust <- c(ave.SilWid.hclust,NA)
        ave.SilWid.PAM <- c(ave.SilWid.PAM,NA)
        ave.SilWid.MC <- c(ave.SilWid.MC,NA)
        cophenetic.corr.hclust <- c(cophenetic.corr.hclust,NA)
        mmcs_pam <- c(mmcs_pam,NA)
        mmcs_hclust <- c(mmcs_hclust,NA)
        mmcs_MC <- c(mmcs_MC,NA)
      }
      
      
      #error
      rmse_list <- c()
      kld_list <- c()
      rmse_orig_list <- c()
      kld_orig_list <- c()
      for (i in 1:saved_nmf_runs){
        selection <- ((i-1)*ns+1):(i*ns)
        current_cat <- boot_cat[[boot_tracker[selection][1]]]
        if(!is.null(matrix_of_fixed_signatures)){
          selection_e <- ((i-1)*(ns+ncol(matrix_of_fixed_signatures))+1):(i*(ns+ncol(matrix_of_fixed_signatures)))
          current_p <- as.matrix(cbind(p_boot[,selection,drop=F],matrix_of_fixed_signatures))
          current_e <- e_boot[,selection_e,drop=F]
          reconstructed_cat <- current_p %*% t(current_e)
        }else{
          current_p <- p_boot[,selection,drop=F]
          current_e <- e_boot[,selection,drop=F]
          reconstructed_cat <- current_p %*% t(current_e)
        }
        rmse_list <- c(rmse_list,sqrt(sum((current_cat - reconstructed_cat)^2)/(dim(current_cat)[1]*dim(current_cat)[2])))
        kld_list <- c(kld_list,KLD(current_cat,reconstructed_cat))
        rmse_orig_list <- c(rmse_orig_list,sqrt(sum((ncat - reconstructed_cat)^2)/(dim(ncat)[1]*dim(ncat)[2])))
        kld_orig_list <- c(kld_orig_list,KLD(ncat,reconstructed_cat))
      }
      # #Save errors
      ave.RMSE <- c(ave.RMSE,mean(rmse_list))
      sd.RMSE <- c(sd.RMSE,sd(rmse_list))
      ave.KLD <- c(ave.KLD,mean(kld_list))
      sd.KLD <- c(sd.KLD,sd(kld_list))
      # 
      ave.RMSE.orig <- c(ave.RMSE.orig,mean(rmse_orig_list))
      sd.RMSE.orig <- c(sd.RMSE.orig,sd(rmse_orig_list))
      ave.KLD.orig <- c(ave.KLD.orig,mean(kld_orig_list))
      sd.KLD.orig <- c(sd.KLD.orig,sd(kld_orig_list))
      
      #use requested medoids
      medoids_final <- medoids_hclust
      if(clusteringMethod=="PAM") medoids_final <- medoids_PAM
      if(clusteringMethod=="MC") medoids_final <- medoids_MC
      
      #plot signatures and compute similarity to known signatures
      clustering_list <- ""
      clustering_list_tag <- clusteringMethod
      if (plotResultsFromAllClusteringMethods) clustering_list <- c("","_HC","_PAM","_MC")
      if (plotResultsFromAllClusteringMethods) clustering_list_tag <- c(clusteringMethod,"HC","PAM","MC")
      
      for (cli in 1:length(clustering_list)){
        cl <- clustering_list[cli]
        cl_tag <- clustering_list_tag[cli]
        signature_names <- paste0("S",1:length(medoids_final))
        
        if (cl==""){
          signature_data_matrix <- p_boot[,medoids_final,drop=FALSE]/matrix(data=rep(apply(p_boot[,medoids_final,drop=FALSE],2,sum),nrow(p_boot)),nrow = nrow(p_boot),byrow = TRUE)
        }else if(cl=="_HC"){
          signature_data_matrix <- p_boot[,medoids_hclust,drop=FALSE]/matrix(data=rep(apply(p_boot[,medoids_hclust,drop=FALSE],2,sum),nrow(p_boot)),nrow = nrow(p_boot),byrow = TRUE)
        }else if(cl=="_PAM"){
          signature_data_matrix <- p_boot[,medoids_PAM,drop=FALSE]/matrix(data=rep(apply(p_boot[,medoids_PAM,drop=FALSE],2,sum),nrow(p_boot)),nrow = nrow(p_boot),byrow = TRUE)
        }else if(cl=="_MC"){
          signature_data_matrix <- p_boot[,medoids_MC,drop=FALSE]/matrix(data=rep(apply(p_boot[,medoids_MC,drop=FALSE],2,sum),nrow(p_boot)),nrow = nrow(p_boot),byrow = TRUE)
        }
        
        colnames(signature_data_matrix) <- signature_names
        row.names(signature_data_matrix) <- row.names(p_boot)
        if(!is.null(matrix_of_fixed_signatures)){
          signature_data_matrix <- cbind(matrix_of_fixed_signatures,signature_data_matrix)
        }
        subs_file <- paste0(outNsDir,"Sigs_plot_",group,"_ns",ns,"_nboots",nboots,cl,".pdf")
        if (type_of_extraction == "subs"){
          #plotSubsSignatures(signature_data_matrix,subs_file,plot_sum = FALSE,overall_title = paste0("Medoids Signatures when extracting ",ns))
          plotSubsSignatures_withMeanSd(signature_data_matrix,mean_signatures[[cl_tag]],sd_signatures[[cl_tag]],subs_file,plot_sum = FALSE,overall_title = paste0("Medoids Signatures when extracting ",ns))
        }else if (type_of_extraction == "rearr"){
          #plotRearrSignatures(signature_data_matrix,subs_file,plot_sum = FALSE,overall_title = paste0("Medoids Signatures when extracting ",ns))
          plotRearrSignatures_withMeanSd(signature_data_matrix,mean_signatures[[cl_tag]],sd_signatures[[cl_tag]],subs_file,plot_sum = FALSE,overall_title = paste0("Medoids Signatures when extracting ",ns))
        }else if (type_of_extraction == "generic"){
          #plotGenericSignatures(signature_data_matrix,subs_file,plot_sum = FALSE,overall_title = paste0("Medoids Signatures when extracting ",ns))
          plotGenericSignatures_withMeanSd(signature_data_matrix,mean_signatures[[cl_tag]],sd_signatures[[cl_tag]],subs_file,plot_sum = FALSE,overall_title = paste0("Medoids Signatures when extracting ",ns))
        }else if (type_of_extraction == "dnv"){
          #plotGenericSignatures(signature_data_matrix,subs_file,plot_sum = FALSE,overall_title = paste0("Medoids Signatures when extracting ",ns))
          plotDNVSignatures_withMeanSd(signature_data_matrix,mean_signatures[[cl_tag]],sd_signatures[[cl_tag]],subs_file,plot_sum = FALSE,overall_title = paste0("Medoids Signatures when extracting ",ns))
        }
        subs_file <- paste0(outNsDir,"Sigs_plot_",group,"_ns",ns,"_nboots",nboots,cl,".tsv")
        write.table(signature_data_matrix,file = subs_file,
                    sep = "\t",quote = FALSE,row.names = TRUE,col.names = TRUE)
        
        #similarity to known signatures
        if(type_of_extraction=="subs"){
          #find the most similar cosmic signatures
          res_cosmic <- findClosestCOSMIC30_withSimilarity(signature_data_matrix) 
          res_cosmicComb <- findClosestCOSMIC30andCombinations_withSimilarity(signature_data_matrix) 
          res_cosmic_table <- data.frame(res_cosmic,res_cosmicComb)
          cosmic_file <- paste0(outNsDir,"Sigs_cosmicSimilar_",group,"_ns",ns,"_nboots",nboots,cl,".tsv")
          write.table(res_cosmic_table,file = cosmic_file,
                      sep = "\t",quote = FALSE,row.names = TRUE,col.names = TRUE)
        }else if (type_of_extraction=="rearr"){
          res_rearrBreast560 <- findClosestRearrSigsBreast560_withSimilarity(signature_data_matrix)
          res_rearrBreast560_table <- data.frame(res_rearrBreast560)
          res_rearrBreast560_file <- paste0(outNsDir,"Sigs_rearrBreast560Similar_",group,"_ns",ns,"_nboots",nboots,cl,".tsv")
          write.table(res_rearrBreast560_table,file = res_rearrBreast560_file,
                      sep = "\t",quote = FALSE,row.names = TRUE,col.names = TRUE)
        }
      }
      
      if (ns>1){
        #More metrics
        #spread of the clusters: minimum within cluster cosine similarity (MinWCCS)
        MinWCCS.hclust <- minWithinClusterCosSim(clustering = cut_res,distMatrix = distMatrix,parallel = TRUE)
        names(MinWCCS.hclust) <- signature_names
        MinWCCS.MC <- minWithinClusterCosSim(clustering = cut_res_MC,distMatrix = distMatrix,parallel = TRUE)
        names(MinWCCS.MC) <- signature_names
        MinWCCS.PAM <- minWithinClusterCosSim(clustering = clustering_result$clustering,distMatrix = distMatrix,parallel = TRUE)
        names(MinWCCS.PAM) <- signature_names
        
        MinWCCS_file <- paste0(outNsDir,"MinWithinClusterCosSim_",group,"_ns",ns,"_nboots",nboots,".tsv")
        write.table(data.frame(MinWCCS.hclust=MinWCCS.hclust,MinWCCS.PAM=MinWCCS.PAM,MinWCCS.MC=MinWCCS.MC),file = MinWCCS_file,
                    sep = "\t",quote = FALSE,row.names = TRUE,col.names = TRUE)
        min.MinWCCS.hclust <- c(min.MinWCCS.hclust,min(MinWCCS.hclust))
        min.MinWCCS.PAM <- c(min.MinWCCS.PAM,min(MinWCCS.PAM))
        min.MinWCCS.MC <- c(min.MinWCCS.MC,min(MinWCCS.MC))
        
        #cluster neighbours: maximum between cluster cosine similarity (MaxBCCS)
        MaxBCCS.hclust <- maxBetweenClustersCosSim(clustering = cut_res,distMatrix = distMatrix,parallel = TRUE)
        colnames(MaxBCCS.hclust) <- signature_names
        row.names(MaxBCCS.hclust) <- signature_names
        MaxBCCS.MC <- maxBetweenClustersCosSim(clustering = cut_res_MC,distMatrix = distMatrix,parallel = TRUE)
        colnames(MaxBCCS.MC) <- signature_names
        row.names(MaxBCCS.MC) <- signature_names
        MaxBCCS.PAM <- maxBetweenClustersCosSim(clustering = clustering_result$clustering,distMatrix = distMatrix,parallel = TRUE)
        colnames(MaxBCCS.PAM) <- signature_names
        row.names(MaxBCCS.PAM) <- signature_names
        
        MaxBCCS.hclust_file <- paste0(outNsDir,"MaxBetweenClusterCosSim.hclust_",group,"_ns",ns,"_nboots",nboots,".tsv")
        write.table(MaxBCCS.hclust,file = MaxBCCS.hclust_file,
                    sep = "\t",quote = FALSE,row.names = TRUE,col.names = TRUE)
        MaxBCCS.PAM_file <- paste0(outNsDir,"MaxBetweenClusterCosSim.PAM_",group,"_ns",ns,"_nboots",nboots,".tsv")
        write.table(MaxBCCS.PAM,file = MaxBCCS.PAM_file,
                    sep = "\t",quote = FALSE,row.names = TRUE,col.names = TRUE)
        MaxBCCS.MC_file <- paste0(outNsDir,"MaxBetweenClusterCosSim.MC_",group,"_ns",ns,"_nboots",nboots,".tsv")
        write.table(MaxBCCS.MC,file = MaxBCCS.MC_file,
                    sep = "\t",quote = FALSE,row.names = TRUE,col.names = TRUE)
        
        max.MaxBCCS.hclust <- c(max.MaxBCCS.hclust,max(MaxBCCS.hclust - diag(nrow(MaxBCCS.hclust)))) 
        max.MaxBCCS.PAM <- c(max.MaxBCCS.PAM,max(MaxBCCS.PAM - diag(nrow(MaxBCCS.PAM))))
        max.MaxBCCS.MC <- c(max.MaxBCCS.MC,max(MaxBCCS.MC - diag(nrow(MaxBCCS.MC)))) 
      }else if (ns==1){
        min.MinWCCS.hclust <- c(min.MinWCCS.hclust,NA)
        min.MinWCCS.PAM <- c(min.MinWCCS.PAM,NA)
        min.MinWCCS.MC <- c(min.MinWCCS.MC,NA)
        max.MaxBCCS.hclust <- c(max.MaxBCCS.hclust,NA) 
        max.MaxBCCS.PAM <- c(max.MaxBCCS.PAM,NA)
        max.MaxBCCS.MC <- c(max.MaxBCCS.MC,NA)
      }      
      #-------------------
      #----- clean this ns
      #-------------------
      rm(p_boot)
      rm(e_boot)
      rm(distMatrix)
      gc()
    }
    
    #Collect metrics
    overall_metrics <- data.frame(nsig,       #1
                                  ave.RMSE,   #2
                                  sd.RMSE,    #3
                                  ave.RMSE.orig, #4
                                  sd.RMSE.orig,  #5
                                  ave.KLD,       #6
                                  sd.KLD,        #7
                                  ave.KLD.orig,  #8
                                  sd.KLD.orig,   #9
                                  ave.SilWid.hclust, #10
                                  ave.SilWid.PAM,    #11
                                  ave.SilWid.MC,     #12
                                  cophenetic.corr.hclust, #13
                                  min.MinWCCS.hclust, #14
                                  min.MinWCCS.PAM,    #15
                                  min.MinWCCS.MC,     #16
                                  max.MaxBCCS.hclust, #17
                                  max.MaxBCCS.PAM,    #18
                                  max.MaxBCCS.MC,     #19
                                  mmcs_hclust,        #20
                                  mmcs_pam,           #21
                                  mmcs_MC)            #22
    
    #write table with all the metrics
    overall_metrics_file <- paste0(outFilePath,"Sigs_OverallMetrics_",group,"_nboots",nboots,".tsv")
    write.table(overall_metrics,file = overall_metrics_file,
                sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE)
    
    whattoplot_hclust <- c(10,14,17,20,13)
    whattoplot_PAM <- c(11,15,18,21)    
    whattoplot_MC <- c(12,16,19,22)
    
    #plot requested medoids
    whattoplot_final <- whattoplot_hclust
    if(clusteringMethod=="PAM") whattoplot_final <- whattoplot_PAM
    if(clusteringMethod=="MC") whattoplot_final <- whattoplot_MC
    
    for (cl in clustering_list){

      overall_metrics_file <- paste0(outFilePath,"Sigs_OverallMetrics_",group,"_nboots",nboots,cl,".jpg")
      #read file for debugging
      # overall_metrics <- read.table(file = overall_metrics_file,sep = "\t",header = TRUE,as.is = TRUE)
      #specify which columns to plot. nsig and ave.RMSE will be used already, just specify others
      # whattoplot <- c(4,6,7,9,11)
      if(cl==""){
        whattoplot <- whattoplot_final
      }else if(cl=="_HC"){
        whattoplot <- whattoplot_hclust
      }else if(cl=="_PAM"){
        whattoplot <- whattoplot_PAM
      }else if(cl=="_MC"){
        whattoplot <- whattoplot_MC
      }
      plotOverallMetrics(overall_metrics,whattoplot,overall_metrics_file,group,nboots,nmfmethod)
    }
    complete <- TRUE
  }
  

  message("\n\nResults can be found in ", outFilePath)
  
  message("\n------------- COMPUTATION  SUCCESSFULLY COMPLETED ------------\n")
  
}

#############################################


