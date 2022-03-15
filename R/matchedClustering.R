
#' Clustering with Matches
#' 
#' Compute the clustering into ns clusters of elements for which a distance matrix (distMatrix) is given, subject to the constraint that groups of ns elements cannot be inthe same cluster. For example, elements from 1 to ns cannot be in the same cluster, and so it is for elements from ns+1 to 2*ns, 2*ns+1 to 3*ns, and so on.
#' 
#' @param distMatrix squared distance matrix
#' @param ns number of clusters
#' @param maxMatch if TRUE, use the maximum match (assignament problem), otherwise use a stable matching (which may not be maximum). Maximum means that the match will have minimal sum of distances
#' @param parallel use multi-threading TRUE/FALSE, uses doMC and foreach packages
#' @param nparallel how many threads at a time
#' @return returns the clusters assigned to each element (there is an element for each column of distMatrix)
#' @keywords clustering constraints stable matching assignment problem
#' @export
#' @examples
#' res.clust <- matchedClustering(distMatrix,10)
matchedClustering <- function(distMatrix,ns,maxMatch=TRUE,parallel=FALSE,nparallel=1){
  nruns <- ncol(distMatrix)/ns
  if (nruns==1) {
    res <- 1:ns
    names(res) <- colnames(distMatrix)
    return(res)
  }
  #the lower the score the better, its a distance
  match_scores <- matrix(NA,ncol = nruns,nrow = nruns)
  match_tables <- list()
  from_coord_to_list <- function(i,j){
    if (i==j) stop("error in from_coord_to_list: i==j not supported in this matrix representation!")
    if(i>j){ #upper triangular representation of symetric matrix, swap indexes!
      tmp <- i
      i <- j
      j <- tmp
    }
    nruns*(nruns-1)/2 - (nruns-i+1)*(nruns-i)/2 + j - i
  }
  # get_from_match_tables <- function(i,j){
  #   match_tables[[from_coord_to_list(i,j)]]
  # }
  if(!parallel){
    for (i in 1:(nruns-1)){
      for (j in (i+1):nruns){
        if (maxMatch){
          match_found <- maxMatchRuns(distMatrix,i,j,ns)
        }else{
          match_found <- matchRuns(distMatrix,i,j,ns)
        }
        match_tables[[from_coord_to_list(i,j)]] <- match_found
        #compute average distance of match
        mean_distance <- 0
        for (s in 1:ns){
          mean_distance <- mean_distance + distMatrix[match_found[s,1],match_found[s,2]]
        }
        mean_distance <- mean_distance/ns
        match_scores[i,j] <- mean_distance
      }
    }
  }else{
    # library(NMF)
    # library(foreach)
    # library(doParallel)
    # library(doMC)
    doParallel::registerDoParallel(nparallel)
    #run initialisation in parallel!
    res_list <- foreach::foreach(i=1:(nruns-1)) %dopar%{
      jlist <- list()
      for (j in (i+1):nruns){
        jlist[[as.character(j)]] <-list()
        if (maxMatch){
          match_found <- maxMatchRuns(distMatrix,i,j,ns)
        }else{
          match_found <- matchRuns(distMatrix,i,j,ns)
        }
        #compute average distance of match
        mean_distance <- 0
        for (s in 1:ns){
          mean_distance <- mean_distance + distMatrix[match_found[s,1],match_found[s,2]]
        }
        mean_distance <- mean_distance/ns
        jlist[[as.character(j)]][["match"]] <- matrix(match_found,ncol = 2) #need to make a copy
        colnames(jlist[[as.character(j)]][["match"]])=colnames(match_found)
        jlist[[as.character(j)]][["mean_distance"]] <- mean_distance
      }
      #return jlist
      jlist
    }
    #reorganise results
    for (i in 1:(nruns-1)){
      for (j in (i+1):nruns){
        match_tables[[from_coord_to_list(i,j)]] <- res_list[[i]][[as.character(j)]][["match"]]
        match_scores[i,j] <- res_list[[i]][[as.character(j)]][["mean_distance"]]
      }
    }
    rm(res_list)
    gc()
  }
  #combine matches
  max_combined <- 2
  combined_matrices <- list()
  combined_matrices_runs <- list()
  check_if_already_found <- function(i,j){
    posi <- -1
    posj <- -1
    p <- 1
    while (p <= length(combined_matrices_runs) & (posi==-1 | posj==-1)){
      if(i %in% combined_matrices_runs[[p]]) posi <- p
      if(j %in% combined_matrices_runs[[p]]) posj <- p
      p <- p + 1
    }
    return(c(posi,posj))
  }
  #get a list of coordinates in order of lowest to highest score
  list_of_best_matches <- lapply(order(t(match_scores),na.last = TRUE)[1:(nruns*(nruns-1)/2)],function(x) c(ceiling(x/nruns),x - (ceiling(x/nruns)-1)*nruns))
  current_best_match_pos <- 1
  if (max_combined==nruns) { #limit case with only two runs
    combined_matrices[[1]] <- match_tables[[from_coord_to_list(1,2)]]
    combined_matrices_runs[[1]] <- list_of_best_matches[[current_best_match_pos]]
  }
  while(max_combined<nruns){
    #get current best match
    current_match <- list_of_best_matches[[current_best_match_pos]]
    i <- current_match[1]
    j <- current_match[2]
    #first check if either i or j have already been found and whether in one or two tables
    res <- check_if_already_found(i,j)
    if(res[1]==-1 & res[2]==-1){
      #neither i nor j have been seen yet, so just add the table to the list
      combined_matrices[[length(combined_matrices)+1]] <- match_tables[[from_coord_to_list(i,j)]]
      combined_matrices_runs[[length(combined_matrices_runs)+1]] <- current_match
      #go to next
      current_best_match_pos <- current_best_match_pos + 1
    }else if(res[1]==-1 & res[2]!=-1){
      #j was found but i not. Merge current match table with pre-existing match table
      combined_matrices[[res[2]]] <- merge(combined_matrices[[res[2]]],match_tables[[from_coord_to_list(i,j)]])
      combined_matrices_runs[[res[2]]] <- c(combined_matrices_runs[[res[2]]],i)
      #update max table size
      newsize <- length(combined_matrices_runs[[res[2]]])
      if(newsize>max_combined) max_combined <- newsize
      #go to next
      current_best_match_pos <- current_best_match_pos + 1
    }else if(res[1]!=-1 & res[2]==-1){
      #i was found but j not. Merge current match table with pre-existing match table
      combined_matrices[[res[1]]] <- merge(combined_matrices[[res[1]]],match_tables[[from_coord_to_list(i,j)]])
      combined_matrices_runs[[res[1]]] <- c(combined_matrices_runs[[res[1]]],j)
      #update max table size
      newsize <- length(combined_matrices_runs[[res[1]]])
      if(newsize>max_combined) max_combined <- newsize
      #go to next
      current_best_match_pos <- current_best_match_pos + 1      
    }else { #that is if(res[1]!=-1 & res[2]!=-1)
      #both runs have been found before
      if(res[1]==res[2]){
        #runs already in the same table, nothing to do, continue
        current_best_match_pos <- current_best_match_pos + 1
      }else{
        #runs are in different tables. Merge current table with one and then the other
        #merge the second table in the list into the first table in the list and remove the second table
        x <- min(res)
        y <- max(res)
        combined_matrices[[x]] <- merge(combined_matrices[[x]],match_tables[[from_coord_to_list(i,j)]])
        combined_matrices[[x]] <- merge(combined_matrices[[x]],combined_matrices[[y]])
        combined_matrices_runs[[x]] <- c(combined_matrices_runs[[x]],combined_matrices_runs[[y]])
        #remove position y
        if(y<length(combined_matrices_runs)){
          #shift left
          for(yy in y:(length(combined_matrices_runs)-1)){
            combined_matrices_runs[[yy]] <- combined_matrices_runs[[yy+1]]
            combined_matrices[[yy]] <- combined_matrices[[yy+1]]
          }
        }
        #set last element to NULL will remove it
        combined_matrices_runs[[length(combined_matrices_runs)]] <- NULL
        combined_matrices[[length(combined_matrices)]] <- NULL
        
        #update max table size
        newsize <- length(combined_matrices_runs[[x]])
        if(newsize>max_combined) max_combined <- newsize
        #go to next
        current_best_match_pos <- current_best_match_pos + 1
      }
    }
  }
  
  #finally, convert to ns partitions/clusters
  clusters <- matrix(rep(1:ns,nruns),nrow = ns)
  final_res <- clusters[order(unlist(combined_matrices[[1]]))]
  names(final_res) <- colnames(distMatrix)
  return(final_res)
}


matchRuns <- function(distMatrix,run1,run2,ns){
  #this function is based on the 
  #Gale-Shapley stable matching algorithm, 1962
  #as a convention, is and js, give the positions in distMatrix
  #i and j give current positions in the is and js vectors
  #and in the match_is and match_js vectors
  is <- 1:ns + (run1-1)*ns
  js <- 1:ns + (run2-1)*ns
  match_is <- rep(NA,ns)
  match_js <- rep(NA,ns)
  #set of i not yet matched
  not_yet_matched <- 1:ns
  #find order of preference in terms of i and j
  preferences_is <- t(apply(distMatrix[is,js],1,order))
  #vector to keep track of current preference tried by each i in is
  current_preference_is <- rep(1,ns)
  while (length(not_yet_matched)>0){
    #first sig not assigned yet
    i <- not_yet_matched[1]
    #first preference not tried yet
    j <- preferences_is[i,current_preference_is[i]]
    #move on to next preference in case i need to try again later
    current_preference_is[i] <- current_preference_is[i] + 1
    if(is.na(match_js[j])){
      #matched!
      match_is[i] <- j
      match_js[j] <- i
      #and remove i from is
      if(length(not_yet_matched)==1){
        #Algorithm stops here with Stable Match!
        not_yet_matched <- c()
      }else{
        not_yet_matched <- not_yet_matched[2:length(not_yet_matched)]
      }
    }else{
      #j is already matched with i1!
      i1 <- match_js[j]
      if(distMatrix[is[i],js[j]] < distMatrix[is[i1],js[j]]){
        #new i is a better match for j than i1
        #free i1, that is add it back to is vector
        not_yet_matched <- c(not_yet_matched,i1)
        #now form new match
        match_is[i] <- j
        match_js[j] <- i
        #and remove i from is
        not_yet_matched <- not_yet_matched[2:length(not_yet_matched)]
        #also free match of i1 from match_is
        #not necessary but cleaner
        match_is[i1] <- NA
      }else{
        #old i1 is better or equal to i as a match for j,
        #so keep as it is, no need to do anything as in
        #the next round i will try a different j
      }
    }
  }
  #return the match
  res <- matrix(data=c(is,js[match_is]),ncol = 2)
  colnames(res) <- c(run1,run2)
  return(res)
}

maxMatchRuns <- function(distMatrix,run1,run2,ns){
  #instance of the assignment problem
  is <- 1:ns + (run1-1)*ns
  js <- 1:ns + (run2-1)*ns
  res_match <- lpSolve::lp.assign(distMatrix[is,js])
  match_is <- apply(res_match$solution,1,which.max)
  #return the match
  res <- matrix(data=c(is,js[match_is]),ncol = 2)
  colnames(res) <- c(run1,run2)
  return(res)
}


