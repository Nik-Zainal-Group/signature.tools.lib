
#------------------------------------------------
#------- Auxiliary functions ---------------------
#------------------------------------------------

#' @importFrom foreach %dopar%

## Generate a random replicate of the cataloge
# This method guarantees the total number of signatures is unchanged
generateRandMuts <- function(x){
  #consider the following method as a replacement
  full_r <- matrix(nrow = dim(x)[1],ncol = dim(x)[2])
  colnames(full_r) <- colnames(x)
  row.names(full_r) <- row.names(x)
  for (i in 1:ncol(x)){
    if(sum(x[,i]>0)){
      samples <- sample(1:nrow(x),size = sum(x[,i]),prob = x[,i]/sum(x[,i]),replace = TRUE)
      r <- unlist(lapply(1:nrow(x),function(p) sum(samples==p)))
    }else{ #no rearrangments found
      r <- x[,i]
    }
    names(r) <- rownames(x)
    full_r[,i] <- r
  }
  return(full_r)
}

## Remove unused Rows and Columns from the Catalogue
## Remove Mutations with small numbers
preprocessCatalgue <- function(d, mut_thr){

  ## Remove Mutations
  nmut <- apply(d, 1, sum)
  nmut <- nmut/sum(nmut)
  pos <- which(nmut<=mut_thr)
  if(length(pos)>0){
    d <- d[-pos,]
  }
  return(d)
}

#' Sort 96-channel Substitution Catalogues
#'
#' This function sorts a matrix of 96-channel Substitution Catalogues,
#' so that the channels (rows) are in the correct order.
#' where each column is a catalogue and each row is a channel.
#' rownames should be set as the name of the channel in the format
#' 5' base[Normal base>Tumour base]3' base, for example A[C>A]A.
#'
#' @param cat catalogues matrix where each column is a catalogue and each row is a channel. Rownames should be set as the name of the channel in the format 5' base[Normal base>Tumour base]3' base, for example A[C>A]A.
#' @return ordered catalogue
#' @export
sortCatalogue <- function(cat){
  all_bp <- c("A", "C", "G", "T")
  pyr_muts <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  nm <- c()
  for(m in pyr_muts){
    for(a in all_bp){
      for(b in all_bp){
        nm <- c(nm, paste(a, "[",m, "]", b, sep=""))
      }
    }
  }
  return(cat[nm,,drop=FALSE])

}



computeCorrelation <- function(x){
  if (ncol(x)==1){
    #if there is only one column, correlation matrix is 1
    return(matrix(1, 1, 1))
  }
  out <- matrix(NA, ncol(x), ncol(x))
  #diagonal is 1
  for(i in 1:ncol(x)){
    out[i,i] <- 1
  }
  colnames(out) <- colnames(x)
  rownames(out) <- colnames(x)
  for(i in 2:ncol(x)){
    for(j in 1:(i-1)){ #up to i-1, diag already set to 1
      #message(i, " ", j)
      out[i,j] <- cos_sim(as.numeric(x[,i]), as.numeric(x[,j]))
      out[j,i] <- out[i,j] #upper triangular is the same
    }
  }
  return(out)
}

#' Compute Correlation (parallel)
#'
#' Compute the correlation of a set of signatures/catalogues according to cosine similarity
#'
#' @param x catalogues/signature matrix where each column is a catalogue/signature and each row is a channel.
#' @param nparallel how many parallel processes to use
#' @param parallel set to TRUE in order to use parallel
#' @return correlation matrix
#' @export
computeCorrelation_parallel <- function(x,nparallel=1,parallel=FALSE){
  if (ncol(x)==1){
    #if there is only one column, correlation matrix is 1
    return(matrix(1, 1, 1))
  }
  out <- matrix(NA, ncol(x), ncol(x))
  #diagonal is 1
  for(i in 1:ncol(x)){
    out[i,i] <- 1
  }
  colnames(out) <- colnames(x)
  rownames(out) <- colnames(x)
  #correlation matrix is symmetric
  if (parallel){
    # library(foreach)
    # library(doParallel)
    # library(doMC)
    doParallel::registerDoParallel(nparallel)
    par_res <- foreach::foreach(i=2:ncol(x)) %dopar% {
      current_res <- c()
      for(j in 1:(i-1)){ #up to i-1, diag already set to 1
        #message(i, " ", j)
        current_res <- c(current_res,cos_sim(as.numeric(x[,i]), as.numeric(x[,j])))
      }
      current_res
    }
    for(i in 2:ncol(x)){
      out[i,1:(i-1)] <- par_res[[i-1]]
      out[1:(i-1),i] <- par_res[[i-1]]
    }
  }else{
    for(i in 2:ncol(x)){
      for(j in 1:(i-1)){ #up to i-1, diag already set to 1
        #message(i, " ", j)
        out[i,j] <- cos_sim(as.numeric(x[,i]), as.numeric(x[,j]))
        out[j,i] <- out[i,j] #upper triangular is the same
      }
    }
  }
  return(out)
}

#' Cosine Similarity
#'
#' Compute the cosine similarity between two vectors, using the formula sum(a*b)/sqrt(sum(a^2)*sum(b^2)).
#'
#' @param a first vector to compare
#' @param b second vector to compare
#' @return cosine similarity
#' @export
cos_sim <- function(a, b){
  return( sum(a*b)/sqrt(sum(a^2)*sum(b^2)) )
}

#function to compute the average within cluster cosine similarity,
#uses as input the distance matrix constructed as 1 - cosSimMatrix
withinClusterCosSim <- function(clustering,distMatrix,parallel){
  clusters <- unique(clustering)
  clusters <- clusters[order(clusters)]
  res <- c()
  if (parallel==TRUE){
    res <- foreach::foreach(i=clusters) %dopar%{
      samplesInCluster <- names(clustering[clustering==i])
      nSamples <- length(samplesInCluster)
      meanCosSim <- 0
      if(nSamples==1){
        meanCosSim <- 1
      }else{
        combinations <- nSamples*(nSamples - 1)/2
        #for efficiency, instead of summing lots of (1 - dist)
        #start from 1*combinations and subtract the dists
        sumCosSim <- combinations
        for(j in 1:(nSamples-1)){
          #for(w in (j+1):(nSamples)){
          sumCosSim <- sumCosSim - sum(distMatrix[samplesInCluster[j],samplesInCluster[1:nSamples > j]])
          #}
        }
        meanCosSim <- sumCosSim/combinations
      }
      meanCosSim
    }
    res <- unlist(res)
  }else{
    for(i in clusters){
      samplesInCluster <- names(clustering[clustering==i])
      nSamples <- length(samplesInCluster)
      combinations <- nSamples*(nSamples - 1)/2
      #for efficiency, instead of summing lots of (1 - dist)
      #start from 1*combinations and subtract the dists
      sumCosSim <- combinations
      for(j in 1:(nSamples-1)){
        #for(w in (j+1):(nSamples)){
        sumCosSim <- sumCosSim - sum(distMatrix[samplesInCluster[j],samplesInCluster[1:nSamples > j]])
        #}
      }
      res <- c(res,sumCosSim/combinations)
    }
  }
  return(res)
}

#function to compute the minimum within cluster cosine similarity,
#uses as input the distance matrix constructed as 1 - cosSimMatrix
minWithinClusterCosSim <- function(clustering,distMatrix,parallel){
  #find max distance and then do minCosSim <- 1 - maxDist
  clusters <- unique(clustering)
  clusters <- clusters[order(clusters)]
  res <- c()
  if (parallel==TRUE){
    res <- foreach::foreach(i=clusters) %dopar%{
      samplesInCluster <- names(clustering[clustering==i])
      nSamples <- length(samplesInCluster)
      maxDist <- 0
      if(nSamples==1){
        #maxDist <- 0
      }else{
        for(j in 1:(nSamples-1)){
          maxDist <- max(maxDist,max(distMatrix[samplesInCluster[j],samplesInCluster[1:nSamples > j]]))
        }
      }
      (1 - maxDist)
    }
    res <- unlist(res)
  }else{
    for(i in clusters){
      samplesInCluster <- names(clustering[clustering==i])
      nSamples <- length(samplesInCluster)
      maxDist <- 0
      if(nSamples==1){
        #maxDist <- 0
      }else{
        for(j in 1:(nSamples-1)){
          maxDist <- max(maxDist,max(distMatrix[samplesInCluster[j],samplesInCluster[1:nSamples > j]]))
        }
      }
      res <- c(res,1 - maxDist)
    }
  }
  return(res)
}

#function to compute the maximum between clusters cosine similarity,
#uses as input the distance matrix constructed as 1 - cosSimMatrix
maxBetweenClustersCosSim <- function(clustering,distMatrix,parallel){
  clusters <- unique(clustering)
  clusters <- clusters[order(clusters)]
  res <- matrix(1,nrow = length(clusters),ncol = length(clusters))
  colnames(res) <- clusters
  rownames(res) <- clusters
  if (parallel==TRUE){
    res_list <- foreach::foreach(c=1:(length(clusters) - 1)) %dopar%{
      i <- clusters[c]
      samplesInCluster <- names(clustering[clustering==i])
      res_table <- data.frame()
      for (d in (c+1):length(clusters)){
        j <- clusters[d]
        samplesInCluster2 <- names(clustering[clustering==j])
        ms <- 1 - min(distMatrix[samplesInCluster,samplesInCluster2])
        res_table <- rbind(res_table,data.frame(r=i,c=j,maxSim=ms))
      }
      res_table
    }

    for (i in 1:length(res_list)){
      currentTable <- res_list[[i]]
      for(j in 1:nrow(currentTable)){
        res[currentTable$r[j],currentTable$c[j]] <- currentTable$maxSim[j]
        res[currentTable$c[j],currentTable$r[j]] <- currentTable$maxSim[j]
      }
    }
  }else{
    for (c in 1:(length(clusters) - 1)){
      i <- clusters[c]
      samplesInCluster <- names(clustering[clustering==i])
      for (d in (c+1):length(clusters)){
        j <- clusters[d]
        samplesInCluster2 <- names(clustering[clustering==j])
        ms <- 1 - min(distMatrix[samplesInCluster,samplesInCluster2])
        res[i,j] <- ms
        res[j,i] <- ms
      }
    }
  }
  return(res)
}

plotHierarchicalCluster <- function(fit_clust,outFilePath,group,ns,nboots){
  output_file <- paste0(outFilePath,"Sigs_Cluster_Average_",group,"_ns",ns,"_nboots",nboots,".jpg")
  jpeg(output_file,width = 12*nboots*ns,height = 500,res = 80)
  plot(fit_clust)
  abline(a=0.1,b=0,col="red")
  dev.off()
}


plotWithinClusterCosSim <- function(cosSimHClust,cosSimPAM,cosSimMC,outFilePath,group,ns,nboots){
  dists <- c(cosSimHClust,cosSimPAM,cosSimMC)
  clustermethod <- c(rep("hclust",length(cosSimHClust)),rep("pam",length(cosSimPAM)),rep("MC",length(cosSimMC)))
  output_file <- paste0(outFilePath,"Sigs_WithinClusterCosSim_",group,"_ns",ns,"_nboots",nboots,".jpg")
  jpeg(output_file,width = 600,height = 500,res = 100)
  boxplot(dists ~ clustermethod, lwd = 2, ylab = 'mean Cosine Similarity',xlab = 'method',
          main = paste0("Within Cluster Cosine Similarity\n",group,", nSig=",ns))
  stripchart( dists ~ clustermethod, vertical = TRUE,
              method = "jitter", add = TRUE, pch = 20, col = 'blue')
  dev.off()
}

plotWithinClusterSilWidth <- function(sil_hclust,sil_pam,sil_MC,outFilePath,group,ns,nboots){
  dists <- c(sil_hclust$clus.avg.widths,sil_pam$clus.avg.widths,sil_MC$clus.avg.widths)
  clustermethod <- c(rep("hclust",length(sil_hclust$clus.avg.widths)),rep("pam",length(sil_pam$clus.avg.widths)),rep("MC",length(sil_pam$clus.avg.widths)))
  output_file <- paste0(outFilePath,"Sigs_WithinClusterSilWidth_",group,"_ns",ns,"_nboots",nboots,".jpg")
  jpeg(output_file,width = 600,height = 500,res = 100)
  boxplot(dists ~ clustermethod, lwd = 2, ylab = 'mean Silhouette Width',xlab = 'method',
          main = paste0("Within Cluster Silhouette Width\n",group,", nSig=",ns))
  stripchart( dists ~ clustermethod, vertical = TRUE,
              method = "jitter", add = TRUE, pch = 20, col = 'blue')
  dev.off()
}


plotOverallMetrics <- function(overall_metrics,whattoplot,overall_metrics_file,group,nboots,nmfmethod){
  jpeg(overall_metrics_file,width = 800,height = 500,res = 100)
  par(mar = c(5, 4, 4, 12),mgp = c(2.5,1,0))
  if(nmfmethod=="lee"){
    max_error <- max(overall_metrics$ave.RMSE,overall_metrics$ave.RMSE.orig)
    plot(overall_metrics$nsig,overall_metrics$ave.RMSE/max_error,type="l",
         ylab = "Score",xlab = "Number of extracted signatures",lwd = 2,
         ylim = c(min(min(overall_metrics$ave.RMSE/max_error,overall_metrics$ave.RMSE.orig/max_error),
                      min(overall_metrics[,whattoplot],na.rm = TRUE))*0.95,1),
         main = paste0("Overall Metrics\n(",group,", bootstraps=",nboots,")"))
    lines(overall_metrics$nsig,overall_metrics$ave.RMSE.orig/max_error,type="l",col="black",lwd = 2,lty = 2)
  }else{
    max_error <- max(overall_metrics$ave.KLD,overall_metrics$ave.KLD.orig)
    plot(overall_metrics$nsig,overall_metrics$ave.KLD/max_error,type="l",
         ylab = "Score",xlab = "Number of extracted signatures",lwd = 2,
         ylim = c(min(min(overall_metrics$ave.KLD/max_error,overall_metrics$ave.KLD.orig/max_error),
                      min(overall_metrics[,whattoplot],na.rm = TRUE))*0.95,1),
         main = paste0("Overall Metrics\n(",group,", bootstraps=",nboots,")"))
    lines(overall_metrics$nsig,overall_metrics$ave.KLD.orig/max_error,type="l",col="black",lwd = 2, lty = 2)
  }
  abline(h = 0.9,lty = 2)
  colours_list <- c("red","green","blue","purple","orange","brown","yellow")
  for(i in 1:length(whattoplot)){
    pos <- whattoplot[i]
    lines(overall_metrics$nsig,overall_metrics[,pos],type="l",col=colours_list[i],lwd = 2)
  }
  # lines(overall_metrics$nsig,overall_metrics$proportion.tooSimilar.Signatures,type="l",col="brown",lwd = 2)
  # legend("right", c("norm.Error","ave.CosSim.hclust","ave.CosSim.PAM","ave.SilWid.hclust","ave.SilWid.PAM","cophenetic.corr.hclust","prop.tooSimilar.Sig"), xpd = TRUE, horiz = FALSE, inset = c(-0.45,-0.4),lty = rep(1,7),
  #        bty = "n", col = c("black","red","green","blue","purple","orange","brown"),lwd = 2, cex = 0.9)
  legend("right", c("norm.Error","norm.Error (orig. cat.)",colnames(overall_metrics)[whattoplot]), xpd = TRUE, horiz = FALSE, inset = c(-0.45,-0.4),lty = c(1,2,rep(1,length(whattoplot))),
         bty = "n", col = c("black","black",colours_list[1:length(whattoplot)]),lwd = 2, cex = 0.9)
  dev.off()
}

findMedoidsHclust <- function(distMatrix,cut_res){
  clusters <- unique(cut_res)
  clusters <- clusters[order(clusters)]
  medoids_hclust <- c()
  for (j in clusters){
    current_cluster <- cut_res[cut_res==j]
    current_cluster_names <- names(current_cluster)
    nInCluster <- length(current_cluster)
    if (nInCluster==1){
      medoids_hclust <- c(medoids_hclust,names(current_cluster))
    }else{
      for (p in 1:nInCluster){
        current_cluster[current_cluster_names[p]] <- sum(distMatrix[current_cluster_names[p],current_cluster_names[-p]])/(nInCluster-1)
      }
      medoids_hclust <- c(medoids_hclust,names(which.min(current_cluster)[1]))
    }
    #message(names(which.min(current_cluster)))
    #message(nInCluster)
    #message(medoids_hclust)
  }
  return(medoids_hclust)
}

average_medoids_cosSim <- function(distMatrix,medoids){
  #at least two clusters!
  nClusters <- length(medoids)
  medoidsNumber <- 1:nClusters
  nCombinations <- nClusters*(nClusters - 1)/2
  amcs <- 0
  for (j in (medoidsNumber-1)){
    amcs <- amcs + sum(distMatrix[medoids[j],medoids[medoidsNumber>j]])
  }
  amcs <- (nCombinations - amcs)/nCombinations
  return(amcs)
}

medoids_cosSimMatrix <- function(p_boot,medoids){
  return(computeCorrelation_parallel(p_boot[,medoids]))
}

computePropTooSimilar <- function(distMatrix,saved_nmf_runs,ns){
  countTooSimilar <- 0
  for (i in 1:length(saved_nmf_runs)){
    if (max(distMatrix[((i-1)*ns + 1):(ns*i),((i-1)*ns + 1):(ns*i)]) > 0.9) countTooSimilar <- countTooSimilar + 1
  }
  return(countTooSimilar/saved_nmf_runs)
}


#' Plot Generic Signatures or Catalogues
#'
#' Function to plot one or more signatures or catalogues with an arbitrary number of channels. Channel names will not be plotted and all channels will be plotted as bars of the same colour.
#'
#' @param signature_data_matrix matrix of signatures, signatures as columns and channels as rows
#' @param output_file set output file, should end with ".jpg" or ".pdf". If output_file==null, output will not be to a file, but will still run the plot functions. The option output_file==null can be used to add this plot to a larger output file.
#' @param plot_sum whether the sum of the channels should be plotted. If plotting signatures this should be FALSE, but if plotting sample catalogues, this can be set to TRUE to display the number of mutations in each sample.
#' @param overall_title set the overall title of the plot
#' @param mar set the option par(mar=mar)
#' @param howManyInOnePage how many signatures or catalogues should be plotted on one page. Multiple pages are plotted if more signatures/catalogues to plot have been requested
#' @param ncolumns how many columns should be used to arrange the signatures/catalogues to plot
#' @export
plotGenericSignatures <- function(signature_data_matrix,
                                  output_file = NULL,
                                  plot_sum = TRUE,
                                  overall_title = "",
                                  add_to_titles = NULL,
                                  mar=NULL,
                                  howManyInOnePage=100,
                                  ncolumns=1){
  if(!is.null(output_file)) plottype <- substr(output_file,nchar(output_file)-2,nchar(output_file))
  npages <- ceiling(ncol(signature_data_matrix)/howManyInOnePage)
  if(!is.null(output_file)) rootoutput_file <- substr(output_file,1,nchar(output_file)-4)
  for(i in 1:npages){
    ifrom <- howManyInOnePage*(i-1) + 1
    ito <- min(ncol(signature_data_matrix),howManyInOnePage*i)
    tmpmatrix <- signature_data_matrix[,ifrom:ito,drop=F]
    if (!is.null(add_to_titles)) tmpadd <- add_to_titles[ifrom:ito]
    if(npages>1 & !is.null(output_file)) output_file <- paste0(rootoutput_file,"_",i,"of",npages,".",plottype)
    # now plot
    nplotrows <- ceiling(ncol(tmpmatrix)/ncolumns)
    if(!is.null(output_file)) {
      if(plottype=="jpg"){
        jpeg(output_file,width = ncolumns*800,height = nplotrows*400,res = 190)
      }else if(plottype=="pdf"){
        pdf(output_file,width = ncolumns*8,height = nplotrows*4+0.5,pointsize = 26)
      }
      par(mfrow = c(nplotrows, ncolumns),omi=c(0,0,0.5,0),cex=0.7)
    }
    for (pos in 1:ncol(tmpmatrix)){
      par(mgp=c(1,1,0))
      if(is.null(mar)){
        par(mar=c(3,3.5,2,1))
      }else{
        par(mar=mar)
      }
      title <- colnames(tmpmatrix)[pos]
      if (!is.null(add_to_titles)) title <- paste0(title," ",tmpadd[pos])
      if (plot_sum) title <- paste0(title,"\n(",round(sum(tmpmatrix[,pos]))," mutations)")

      barplot(tmpmatrix[,pos],
              main = title,
              names.arg = c(rep("",nrow(tmpmatrix))),
              col="skyblue",
              beside = TRUE,
              xlab = "channels",
              cex.main = 0.9,
              cex.names = 1)
    }
    title(main = overall_title,outer = TRUE,cex.main = 1.5)
    if(!is.null(output_file)) dev.off()
  }
}

plotGenericSignatures_withMeanSd <- function(signature_data_matrix,
                                             mean_matrix,
                                             sd_matrix,
                                             output_file = NULL,
                                             plot_sum = TRUE,
                                             overall_title = "",
                                             add_to_titles = NULL,
                                             mar=NULL){
  if(!is.null(output_file)) plottype <- substr(output_file,nchar(output_file)-2,nchar(output_file))
  # rearr.colours <- c(rep("blue",16),rep("black",16),rep("red",16),rep("grey",16),rep("green",16),rep("pink",16))
  nplotrows <- ncol(signature_data_matrix)
  if(!is.null(output_file)) {
    if(plottype=="jpg"){
      jpeg(output_file,width = 2*800,height = nplotrows*400,res = 190)
    }else if (plottype=="pdf"){
      pdf(output_file,width = 2*8,height = nplotrows*4+0.5,pointsize = 26)
    }
  }
  par(mfrow = c(nplotrows, 2),omi=c(0,0,0.5,0))
  par(mgp=c(1,1,0))
  if(is.null(mar)){
    par(mar=c(3,3,2,1))
  }else{
    par(mar=mar)
  }
  for (pos in 1:ncol(signature_data_matrix)){
    ylimit <- c(0,max(signature_data_matrix[,pos],mean_matrix[,pos]+sd_matrix[,pos]))
    title <- colnames(signature_data_matrix)[pos]
    if (!is.null(add_to_titles)) title <- paste0(title," ",add_to_titles[pos])
    if (plot_sum) title <- paste0(title,"\n(",round(sum(signature_data_matrix[,pos]))," mutations)")
    barplot(signature_data_matrix[,pos],
            main = title,
            names.arg = c(rep("",nrow(signature_data_matrix))),
            col="skyblue",
            beside = TRUE,
            xlab = "channels",
            cex.main = 0.9,
            ylim = ylimit,
            cex.names = 1)
    barCenters <- barplot(mean_matrix[,pos],
                          main = "mean and sd of cluster",
                          names.arg = c(rep("",nrow(signature_data_matrix))),
                          col="skyblue",
                          beside = TRUE,
                          xlab = "channels",
                          cex.main = 0.9,
                          ylim = ylimit,
                          cex.names = 1)
    segments(barCenters, mean_matrix[,pos] - sd_matrix[,pos], barCenters,
             mean_matrix[,pos] + sd_matrix[,pos], lwd = 1.5)

    arrows(barCenters, mean_matrix[,pos] - sd_matrix[,pos], barCenters,
           mean_matrix[,pos] + sd_matrix[,pos], lwd = 1.5, angle = 90,
           code = 3, length = 0.05)
  }
  title(main = overall_title,outer = TRUE,cex.main = 1.5)
  if(!is.null(output_file)) dev.off()
}

#' Plot Substitution Signatures or Catalogues
#'
#' Function to plot one or more substitution signatures or catalogues.
#'
#' @param signature_data_matrix matrix of signatures, signatures as columns and channels as rows
#' @param output_file set output file, should end with ".jpg" or ".pdf". If output_file==null, output will not be to a file, but will still run the plot functions. The option output_file==null can be used to add this plot to a larger output file.
#' @param plot_sum whether the sum of the channels should be plotted. If plotting signatures this should be FALSE, but if plotting sample catalogues, this can be set to TRUE to display the number of mutations in each sample.
#' @param overall_title set the overall title of the plot
#' @param add_to_titles text to be added to the titles of each catalogue plot
#' @param mar set the option par(mar=mar)
#' @param howManyInOnePage how many signatures or catalogues should be plotted on one page. Multiple pages are plotted if more signatures/catalogues to plot have been requested
#' @param ncolumns how many columns should be used to arrange the signatures/catalogues to plot
#' @export
plotSubsSignatures <- function(signature_data_matrix,
                               output_file = NULL,
                               plot_sum = TRUE,
                               overall_title = "",
                               add_to_titles = NULL,
                               mar=NULL,
                               howManyInOnePage=100,
                               ncolumns=1){
  # colnames(signature_data_matrix) <- sapply(colnames(signature_data_matrix),function(x) if (nchar(x)>30) paste0(substr(x,1,23),"...") else x)
  # plotcolours <- c("blue","black","red","gray","green","pink")
  plotcolours <- c(rgb(5,195,239,maxColorValue = 255),
                   rgb(0,0,0,maxColorValue = 255),
                   rgb(230,47,41,maxColorValue = 255),
                   rgb(208,207,207,maxColorValue = 255),
                   rgb(169,212,108,maxColorValue = 255),
                   rgb(238,205,204,maxColorValue = 255))
  if(!is.null(output_file)) plottype <- substr(output_file,nchar(output_file)-2,nchar(output_file))
  rearr.colours <- c(rep(plotcolours[1],16),rep(plotcolours[2],16),rep(plotcolours[3],16),rep(plotcolours[4],16),rep(plotcolours[5],16),rep(plotcolours[6],16))
  npages <- ceiling(ncol(signature_data_matrix)/howManyInOnePage)
  if(!is.null(output_file)) rootoutput_file <- substr(output_file,1,nchar(output_file)-4)
  for(i in 1:npages){
    ifrom <- howManyInOnePage*(i-1) + 1
    ito <- min(ncol(signature_data_matrix),howManyInOnePage*i)
    tmpmatrix <- signature_data_matrix[,ifrom:ito,drop=F]
    if (!is.null(add_to_titles)) tmpadd <- add_to_titles[ifrom:ito]
    if(npages>1 & !is.null(output_file)) output_file <- paste0(rootoutput_file,"_",i,"of",npages,".",plottype)
    # now plot
    nplotrows <- ceiling(ncol(tmpmatrix)/ncolumns)
    if(!is.null(output_file)) {
      if(plottype=="jpg"){
        jpeg(output_file,width = ncolumns*800,height = nplotrows*300,res = 220)
      }else if(plottype=="pdf"){
        pdf(output_file,width = ncolumns*8,height = nplotrows*3+0.5,pointsize = 26)
      }
      par(mfrow = c(nplotrows, ncolumns),omi=c(0,0,0.5,0),cex=0.7)
    }
    for (pos in 1:ncol(tmpmatrix)){
      if(is.null(mar)){
        par(mar=c(2,3.5,2,1))
      }else{
        par(mar=mar)
      }
      title <- colnames(tmpmatrix)[pos]
      if (!is.null(add_to_titles)) title <- paste0(title," ",tmpadd[pos])
      if (plot_sum) title <- paste0(title,"\n(",round(sum(tmpmatrix[,pos]))," SNVs)")
      muttypes <- c("C>A","C>G","C>T","T>A","T>C","T>G")
      xlabels <- rep("",96)
      barplot(tmpmatrix[,pos],
              main = title,
              names.arg = xlabels,
              col=rearr.colours,
              beside = TRUE,
              las=2,
              cex.main = 0.9,
              cex.names = 1,border = NA,space = 0.2)
      par(xpd=TRUE)
      par(usr = c(0, 1, 0, 1))
      recttop <- -0.02
      rectbottom <- -0.16
      start1 <- 0.035
      gap <- 0.155
      rect(start1, rectbottom, start1+gap, recttop,col = plotcolours[1],border = NA)
      rect(start1+gap, rectbottom, start1+2*gap, recttop,col = plotcolours[2],border = NA)
      rect(start1+2*gap, rectbottom, start1+3*gap, recttop,col = plotcolours[3],border = NA)
      rect(start1+3*gap, rectbottom, start1+4*gap, recttop,col = plotcolours[4],border = NA)
      rect(start1+4*gap, rectbottom, start1+5*gap, recttop,col = plotcolours[5],border = NA)
      rect(start1+5*gap, rectbottom, start1+6*gap, recttop,col = plotcolours[6],border = NA)
      textposx <- 0.04+seq(8,88,16)/104
      text(x = textposx[1:3],y = -0.09,labels = muttypes[1:3],col = "white",font = 2)
      text(x = textposx[4:6],y = -0.09,labels = muttypes[4:6],col = "black",font = 2)
      #shadowtext(x = 0.04+seq(8,88,16)/104,y = rep(-0.09,6),labels = muttypes,col = "white",bg = "black",r=0.2)
      par(xpd=FALSE)
    }
    title(main = overall_title,outer = TRUE,cex.main = 1.5)
    if(!is.null(output_file)) dev.off()
  }
}

plotSubsSignatures_withMeanSd <- function(signature_data_matrix,
                                          mean_matrix,
                                          sd_matrix,
                                          output_file = NULL,
                                          plot_sum = TRUE,
                                          overall_title = "",
                                          add_to_titles = NULL,
                                          mar=NULL){
  # colnames(signature_data_matrix) <- sapply(colnames(signature_data_matrix),function(x) if (nchar(x)>30) paste0(substr(x,1,22),"...") else x)
  # plotcolours <- c("blue","black","red","gray","green","pink")
  plotcolours <- c(rgb(5,195,239,maxColorValue = 255),
                   rgb(0,0,0,maxColorValue = 255),
                   rgb(230,47,41,maxColorValue = 255),
                   rgb(208,207,207,maxColorValue = 255),
                   rgb(169,212,108,maxColorValue = 255),
                   rgb(238,205,204,maxColorValue = 255))
  if(!is.null(output_file)) plottype <- substr(output_file,nchar(output_file)-2,nchar(output_file))
  rearr.colours <- c(rep(plotcolours[1],16),rep(plotcolours[2],16),rep(plotcolours[3],16),rep(plotcolours[4],16),rep(plotcolours[5],16),rep(plotcolours[6],16))
  nplotrows <- ncol(signature_data_matrix)
  if(!is.null(output_file)) {
    if(plottype=="jpg"){
      jpeg(output_file,width = 2*800,height = nplotrows*300,res = 220)
    }else if(plottype=="pdf"){
      pdf(output_file,width = 2*8,height = nplotrows*3+0.5,pointsize = 26)
    }
  }
  par(mfrow = c(nplotrows, 2),omi=c(0,0,0.5,0),cex=0.7)
  if(is.null(mar)){
    par(mar=c(2,3,2,1))
  }else{
    par(mar=mar)
  }
  muttypes <- c("C>A","C>G","C>T","T>A","T>C","T>G")
  xlabels <- rep("",96)
  # xlabels[8] <- "C > A"
  # xlabels[24] <- "C > G"
  # xlabels[40] <- "C > T"
  # xlabels[56] <- "T > A"
  # xlabels[72] <- "T > C"
  # xlabels[88] <- "T > G"
  for (pos in 1:ncol(signature_data_matrix)){
    ylimit <- c(0,max(signature_data_matrix[,pos],mean_matrix[,pos]+sd_matrix[,pos]))
    title <- colnames(signature_data_matrix)[pos]
    if (!is.null(add_to_titles)) title <- paste0(title," ",add_to_titles[pos])
    if (plot_sum) title <- paste0(title,"\n(",round(sum(signature_data_matrix[,pos]))," SNVs)")
    barplot(signature_data_matrix[,pos],
            main = title,
            #names.arg = row.names(signature_data_matrix),
            names.arg = xlabels,
            col=rearr.colours,
            beside = TRUE,
            ylim = ylimit,
            las=2,
            cex.main = 0.9,
            cex.names = 1)
    par(xpd=TRUE)
    par(usr = c(0, 1, 0, 1))
    recttop <- -0.02
    rectbottom <- -0.16
    start1 <- 0.035
    gap <- 0.155
    rect(start1, rectbottom, start1+gap, recttop,col = plotcolours[1],border = NA)
    rect(start1+gap, rectbottom, start1+2*gap, recttop,col = plotcolours[2],border = NA)
    rect(start1+2*gap, rectbottom, start1+3*gap, recttop,col = plotcolours[3],border = NA)
    rect(start1+3*gap, rectbottom, start1+4*gap, recttop,col = plotcolours[4],border = NA)
    rect(start1+4*gap, rectbottom, start1+5*gap, recttop,col = plotcolours[5],border = NA)
    rect(start1+5*gap, rectbottom, start1+6*gap, recttop,col = plotcolours[6],border = NA)
    textposx <- 0.04+seq(8,88,16)/104
    text(x = textposx[1:3],y = -0.09,labels = muttypes[1:3],col = "white",font = 2)
    text(x = textposx[4:6],y = -0.09,labels = muttypes[4:6],col = "black",font = 2)
    par(xpd=FALSE)
    barCenters <- barplot(mean_matrix[,pos],
                          main = "mean and sd of cluster",
                          #names.arg = row.names(signature_data_matrix),
                          names.arg = xlabels,
                          col=rearr.colours,
                          #border = NA,
                          beside = TRUE,
                          ylim = ylimit,
                          las=2,
                          cex.main = 0.9,
                          cex.names = 1,border = NA,space = 0.2)
    # segments(barCenters, mean_matrix[,pos] - sd_matrix[,pos], barCenters,
    #          mean_matrix[,pos] + sd_matrix[,pos], lwd = 1.5)
    segments(barCenters, mean_matrix[,pos], barCenters,
             mean_matrix[,pos] + sd_matrix[,pos], lwd = 1)
    # arrows(barCenters, mean_matrix[,pos] - sd_matrix[,pos], barCenters,
    #        mean_matrix[,pos] + sd_matrix[,pos], lwd = 1.5, angle = 90,
    #        code = 3, length = 0.05)
    par(xpd=TRUE)
    par(usr = c(0, 1, 0, 1))
    recttop <- -0.02
    rectbottom <- -0.16
    start1 <- 0.035
    gap <- 0.155
    rect(start1, rectbottom, start1+gap, recttop,col = plotcolours[1],border = NA)
    rect(start1+gap, rectbottom, start1+2*gap, recttop,col = plotcolours[2],border = NA)
    rect(start1+2*gap, rectbottom, start1+3*gap, recttop,col = plotcolours[3],border = NA)
    rect(start1+3*gap, rectbottom, start1+4*gap, recttop,col = plotcolours[4],border = NA)
    rect(start1+4*gap, rectbottom, start1+5*gap, recttop,col = plotcolours[5],border = NA)
    rect(start1+5*gap, rectbottom, start1+6*gap, recttop,col = plotcolours[6],border = NA)
    textposx <- 0.04+seq(8,88,16)/104
    text(x = textposx[1:3],y = -0.09,labels = muttypes[1:3],col = "white",font = 2)
    text(x = textposx[4:6],y = -0.09,labels = muttypes[4:6],col = "black",font = 2)
    par(xpd=FALSE)
  }
  title(main = overall_title,outer = TRUE,cex.main = 1.5)
  if(!is.null(output_file)) dev.off()
}

#' Plot Rearrangement Signatures or Catalogues
#'
#' Function to plot one or more rearrangement signatures or catalogues.
#'
#' @param signature_data_matrix matrix of signatures, signatures as columns and channels as rows
#' @param output_file set output file, should end with ".jpg" or ".pdf". If output_file==null, output will not be to a file, but will still run the plot functions. The option output_file==null can be used to add this plot to a larger output file.
#' @param plot_sum whether the sum of the channels should be plotted. If plotting signatures this should be FALSE, but if plotting sample catalogues, this can be set to TRUE to display the number of mutations in each sample.
#' @param overall_title set the overall title of the plot
#' @param mar set the option par(mar=mar)
#' @param howManyInOnePage how many signatures or catalogues should be plotted on one page. Multiple pages are plotted if more signatures/catalogues to plot have been requested
#' @param ncolumns how many columns should be used to arrange the signatures/catalogues to plot
#' @export
plotRearrSignatures <-function(signature_data_matrix,
                               output_file = NULL,
                               plot_sum = TRUE,
                               overall_title = "",
                               add_to_titles = NULL,
                               mar=NULL,
                               howManyInOnePage=100,
                               ncolumns=1){
  #This function plots a set of signatures in a single file, three signatures for each row.
  #signature_data_matrix is a data frame that contains the rearrangement signatures.
  #                      The columns are the signatures, while the rows are the 32 features
  # colnames(signature_data_matrix) <- sapply(colnames(signature_data_matrix),function(x) if (nchar(x)>30) paste0(substr(x,1,22),"...") else x)
  if(!is.null(output_file)) plottype <- substr(output_file,nchar(output_file)-2,nchar(output_file))
  del_col = rgb(228,26,28, maxColorValue = 255)
  td_col = rgb(77,175,74, maxColorValue =255)
  inv_col  = rgb(55,126,184, maxColorValue = 255)
  transloc_col = rgb(152,78,163, maxColorValue =255)
  non_clust_col = rgb(240,240,240, maxColorValue =255)
  #rearr.colours <- c(rep("darkblue",16),rep("red",16))
  rearr.colours <- rep(c(rep(del_col,5),rep(td_col,5),rep(inv_col,5),transloc_col),2)
  npages <- ceiling(ncol(signature_data_matrix)/howManyInOnePage)
  if(!is.null(output_file)) rootoutput_file <- substr(output_file,1,nchar(output_file)-4)

  for(i in 1:npages){
    ifrom <- howManyInOnePage*(i-1) + 1
    ito <- min(ncol(signature_data_matrix),howManyInOnePage*i)
    tmpmatrix <- signature_data_matrix[,ifrom:ito,drop=F]
    if (!is.null(add_to_titles)) tmpadd <- add_to_titles[ifrom:ito]
    if(npages>1 & !is.null(output_file)) output_file <- paste0(rootoutput_file,"_",i,"of",npages,".",plottype)
    # now plot
    nplotrows <- ceiling(ncol(tmpmatrix)/ncolumns)
    if(!is.null(output_file)){
      if(plottype=="jpg"){
        jpeg(output_file,width = ncolumns*800,height = nplotrows*500,res = 220)
      }else if(plottype=="pdf"){
        pdf(output_file,width = ncolumns*8,height = nplotrows*5+0.5,pointsize = 26)
      }
      par(mfrow = c(nplotrows, ncolumns),omi=c(0,0,0.5,0),cex=0.7)
    }
    sizes <- c("1-10Kb",
               "10-100Kb",
               "100Kb-1Mb",
               "1Mb-10Mb",
               ">10Mb")
    sizes_names <- c(rep(sizes,3),"",rep(sizes,3),"")
    for (pos in 1:ncol(tmpmatrix)){
      if(is.null(mar)){
        par(mar=c(8,3,2,1))
      }else{
        par(mar=mar)
      }
      title <- colnames(tmpmatrix)[pos]
      if (!is.null(add_to_titles)) title <- paste0(title," ",tmpadd[pos])
      if (plot_sum) title <- paste0(title,"\n(",sum(tmpmatrix[,pos])," SVs)")
      pos <- barplot(tmpmatrix[,pos],
                     main = title,
                     names.arg = NA,
                     #names.arg = sizes_names,
                     col=rearr.colours,
                     beside = FALSE,
                     #las=2,
                     cex.names = 0.8,
                     cex.main = 0.9,
                     #mgp=c(3,2,0),
                     border = 0,
                     space = 0.1)
      axis(1,
           las=2,
           #hadj=0.5,
           at=pos,
           lab=sizes_names,
           #mgp=c(3,2,0),
           col = "transparent",
           line = 1,
           cex.axis = 0.8)
      #save old plot coordinates
      op <- par("usr")
      #set new coordinates
      par(usr = c(0, 1, 0, 1))
      #add graphics
      par(xpd=TRUE)
      start1 <- 0.035
      xsep = 0.145
      start1_text <- 0.11
      tr_size <- 0.03
      #rect(0.1, 0.1, 0.2, 0.2,col = "blue",lwd = 0)
      #rect(0.1, 0.1, 0.11, 0.11,col = "red",lwd = 0)
      stop <- start1
      for(i in 1:2){
        start <- stop
        stop <- start + xsep
        rect(start, -0.14, stop, -0.02,col = del_col,lwd = 0,border = NA)
        text(x = start+0.5*xsep,y = -0.08,"del",col = "white")
        start <- stop
        stop <- start + xsep
        rect(start, -0.14, stop, -0.02,col = td_col,lwd = 0,border = NA)
        text(x = start+0.5*xsep,y = -0.08,"tds",col = "white")
        start <- stop
        stop <- start + xsep
        rect(start, -0.14, stop, -0.02,col = inv_col,lwd = 0,border = NA)
        text(x = start+0.5*xsep,y = -0.08,"inv",col = "white")
        start <- stop
        stop <- start + tr_size
        rect(start, -0.14, stop, -0.02,col = transloc_col,lwd = 0,border = NA)
        text(x = start+0.5*tr_size,y = -0.08,"tr",col = "white")
      }
      xsep2 <- 3*xsep+tr_size
      rect(start1, -0.26, start1+xsep2, -0.14,col = "black",lwd = 0,border = NA)
      text(x = start1+0.5*xsep2,y = -0.2,"clustered",col = "white")
      rect(start1+xsep2, -0.26, start1+2*xsep2, -0.14,col = non_clust_col,lwd = 0,border = NA)
      text(x = start1+1.5*xsep2,y = -0.2,"non-clustered",col = "black")

      #restore old coordinates
      par(usr = op)
    }
    title(main = overall_title,outer = TRUE,cex.main = 1.5)
    if(!is.null(output_file)) dev.off()
  }
}


plotRearrSignatures_withMeanSd <-function(signature_data_matrix,
                                          mean_matrix,
                                          sd_matrix,
                                          output_file = NULL,
                                          plot_sum = TRUE,
                                          overall_title = "",
                                          add_to_titles = NULL,
                                          mar=NULL){
  #This function plots a set of signatures in a single file, three signatures for each row.
  #signature_data_matrix is a data frame that contains the rearrangement signatures.
  #                      The columns are the signatures, while the rows are the 32 features
  # colnames(signature_data_matrix) <- sapply(colnames(signature_data_matrix),function(x) if (nchar(x)>30) paste0(substr(x,1,22),"...") else x)
  if(!is.null(output_file)) plottype <- substr(output_file,nchar(output_file)-2,nchar(output_file))
  del_col = rgb(228,26,28, maxColorValue = 255)
  td_col = rgb(77,175,74, maxColorValue =255)
  inv_col  = rgb(55,126,184, maxColorValue = 255)
  transloc_col = rgb(152,78,163, maxColorValue =255)
  non_clust_col = rgb(240,240,240, maxColorValue =255)
  #rearr.colours <- c(rep("darkblue",16),rep("red",16))
  rearr.colours <- rep(c(rep(del_col,5),rep(td_col,5),rep(inv_col,5),transloc_col),2)
  nplotrows <- ncol(signature_data_matrix)
  if(!is.null(output_file)) {
    if(plottype=="jpg"){
      jpeg(output_file,width = 2*800,height = nplotrows*500,res = 220)
    }else if (plottype=="pdf"){
      pdf(output_file,width = 2*8,height = nplotrows*5+0.5,pointsize = 26)
    }
  }
  par(mfrow = c(nplotrows, 2),omi=c(0,0,0.5,0),cex=0.7)
  sizes <- c("1-10Kb",
             "10-100Kb",
             "100Kb-1Mb",
             "1Mb-10Mb",
             ">10Mb")
  sizes_names <- c(rep(sizes,3),"",rep(sizes,3),"")

  rearrAxis <- function(barCenters,sizes_names){
    axis(1,
         las=2,
         #hadj=0.5,
         at=barCenters,
         lab=sizes_names,
         #mgp=c(3,2,0),
         col = "transparent",
         line = 1,
         cex.axis = 0.8)
    #save old plot coordinates
    op <- par("usr")
    #set new coordinates
    par(usr = c(0, 1, 0, 1))
    #add graphics
    par(xpd=TRUE)
    start1 <- 0.035
    xsep = 0.145
    start1_text <- 0.11
    tr_size <- 0.03
    #rect(0.1, 0.1, 0.2, 0.2,col = "blue",lwd = 0)
    #rect(0.1, 0.1, 0.11, 0.11,col = "red",lwd = 0)
    stop <- start1
    for(i in 1:2){
      start <- stop
      stop <- start + xsep
      rect(start, -0.14, stop, -0.02,col = del_col,lwd = 0,border = NA)
      text(x = start+0.5*xsep,y = -0.08,"del",col = "white")
      start <- stop
      stop <- start + xsep
      rect(start, -0.14, stop, -0.02,col = td_col,lwd = 0,border = NA)
      text(x = start+0.5*xsep,y = -0.08,"tds",col = "white")
      start <- stop
      stop <- start + xsep
      rect(start, -0.14, stop, -0.02,col = inv_col,lwd = 0,border = NA)
      text(x = start+0.5*xsep,y = -0.08,"inv",col = "white")
      start <- stop
      stop <- start + tr_size
      rect(start, -0.14, stop, -0.02,col = transloc_col,lwd = 0,border = NA)
      text(x = start+0.5*tr_size,y = -0.08,"tr",col = "white")
    }
    xsep2 <- 3*xsep+tr_size
    rect(start1, -0.26, start1+xsep2, -0.14,col = "black",lwd = 0,border = NA)
    text(x = start1+0.5*xsep2,y = -0.2,"clustered",col = "white")
    rect(start1+xsep2, -0.26, start1+2*xsep2, -0.14,col = non_clust_col,lwd = 0,border = NA)
    text(x = start1+1.5*xsep2,y = -0.2,"non-clustered",col = "black")

    #restore old coordinates
    par(usr = op)
  }

  for (pos in 1:ncol(signature_data_matrix)){
    ylimit <- c(0,max(signature_data_matrix[,pos],mean_matrix[,pos]+sd_matrix[,pos]))
    if(is.null(mar)){
      par(mar=c(8,3,2,1))
    }else{
      par(mar=mar)
    }
    title <- colnames(signature_data_matrix)[pos]
    if (!is.null(add_to_titles)) title <- paste0(title," ",add_to_titles[pos])
    if (plot_sum) title <- paste0(title,"\n(",sum(signature_data_matrix[,pos])," SVs)")
    barCenters <- barplot(signature_data_matrix[,pos],
                          main = title,
                          names.arg = NA,
                          #names.arg = sizes_names,
                          col=rearr.colours,
                          beside = FALSE,
                          #las=2,
                          cex.names = 0.8,
                          #mgp=c(3,2,0),
                          border = 0,
                          ylim = ylimit,
                          cex.main = 0.9,
                          space = 0.1)
    rearrAxis(barCenters,sizes_names)
    barCenters <- barplot(mean_matrix[,pos],
                          main = "mean and sd of cluster",
                          #names.arg = row.names(signature_data_matrix),
                          names.arg = NA,
                          col=rearr.colours,
                          #border = NA,
                          beside = FALSE,
                          ylim = ylimit,
                          cex.main = 0.9,
                          las=2,
                          border = 0,
                          space = 0.1)
    # segments(barCenters, mean_matrix[,pos] - sd_matrix[,pos], barCenters,
    #          mean_matrix[,pos] + sd_matrix[,pos], lwd = 1.5)
    segments(barCenters, mean_matrix[,pos], barCenters,
             mean_matrix[,pos] + sd_matrix[,pos], lwd = 1.5)
    # arrows(barCenters, mean_matrix[,pos] - sd_matrix[,pos], barCenters,
    #        mean_matrix[,pos] + sd_matrix[,pos], lwd = 1.5, angle = 90,
    #        code = 3, length = 0.05)
    rearrAxis(barCenters,sizes_names)
  }
  title(main = overall_title,outer = TRUE,cex.main = 1.5)
  if(!is.null(output_file)) dev.off()
}


plotCosSimMatrix <- function(CosSimMatrix,output_file,dpi=300,xlabel = "",ylabel = "",thresholdMark = 0.9,extraWidth = 500,extraHeight = 500,ndigitsafterzero = 2){
  # library("ggplot2")

  # Set up the vectors
  signatures.names <- colnames(CosSimMatrix)
  sample.names <- row.names(CosSimMatrix)

  # Create the data frame
  df <- expand.grid(sample.names,signatures.names)
  df$value <- unlist(CosSimMatrix)
  df$labels <- sprintf(paste0("%.",ndigitsafterzero,"f"), df$value)
  df$labels[df$value==0] <- ""

  #Plot the Data (500+150*nsamples)x1200
  g <- ggplot2::ggplot(df, ggplot2::aes(Var1, Var2)) + ggplot2::geom_point(ggplot2::aes(size = value, colour = value>thresholdMark)) + ggplot2::theme_bw() + ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel)
  g <- g + ggplot2::scale_size_continuous(range=c(0,10)) + ggplot2::geom_text(ggplot2::aes(label = labels))
  g + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size=14),
            axis.text.y = ggplot2::element_text(vjust = 1, size=14)) + ggplot2::theme(legend.position="none")
  w <- (extraWidth+150*length(sample.names))/dpi
  h <- (extraHeight+150*length(signatures.names))/dpi
  ggplot2::ggsave(filename = output_file,dpi = dpi,height = h,width = w,limitsize = FALSE)
}

#' plotCosSimSignatures
#'
#' Plot a matrix of cosine similarities between two sets of signatures.
#'
#' @param sig1 matrix with signatures as columns
#' @param sig2 matrix with signatures as columns
#' @param output_file name of the output file (pdf), optional
#' @param cex.numbers scale the text used for the numbers in the matrix
#' @param circlesColBasic colour used for the circles
#' @param circlesColHighlight colour used for the circles that pass the thresholdMark
#' @export
plotCosSimSignatures <- function(sig1,sig2,output_file = NULL,
                                 thresholdMark = NULL,
                                 cex.numbers = 0.7,
                                 circlesColBasic = "#A1CAF1",
                                 circlesColHighlight = "#F6A600"){
  sigcor <- computeCorrelationOfTwoSetsOfSigs(sig1,sig2)
  plotMatrix(sigcor,output_file = output_file,
             thresholdMark = thresholdMark,
             cex.numbers = cex.numbers,
             circlesColBasic = circlesColBasic,
             circlesColHighlight = circlesColHighlight)
}

#' Find closest COSMIC30 signatures
#'
#' Compares a set of signatures to the COSMIC30 signatures and
#' returns the list of signatures identified. For example,
#' c("C1","C3","N1","C13","N2") means that Cosmic (C) signatures 1, 3 and 13 were found, while
#' signatures N1 and N2 are unknown signatures (N for not found), based on a similarity threshold (similarity >=threshold)
#' comparing to COSMIC30
#'
#' @param sigs matrix with 96-channel substitution signatures as columns
#' @param threshold cosine similarity threshold
#' @return the list of signatures identified
#' @export
findClosestCOSMIC30 <- function(sigs,threshold){
  #load COSMIC30
  # cosmic30 <- read.table("../data/COSMIC_signatures.txt", sep="\t", header=T, as.is=T, check.names = FALSE)
  #compute cos sim matrix
  cos_sim_df <- data.frame()
  for (s in colnames(sigs)){
    for(a in colnames(cosmic30)){
      cos_sim_df[s,a] <- cos_sim(sigs[,s],cosmic30[,a])
    }
  }
  max.sim <- apply(cos_sim_df,1,max)
  closestCosmic <- apply(cos_sim_df,1,which.max)
  notFound <- max.sim < threshold
  closestCosmic[notFound] <- paste0("N",1:sum(notFound))
  closestCosmic[!notFound] <- paste0("C",closestCosmic[!notFound])
  return(closestCosmic)
}

#automatically detect similarity with sum of two COSMIC30

#' Find closest COSMIC30 signatures or combination of COSMIC30
#'
#' Compares a set of signatures to the COSMIC30 signatures or the simple sum of all combinations of two COSMIC signatures and
#' returns the list of signatures identified. For example,
#' c("C1","C3","N1","C2+13","N2") means that Cosmic (C) signatures 1, 3 and 2+13 were found, while
#' signatures N1 and N2 are unknown signatures (N for not found), based on a similarity threshold (similarity >=threshold)
#' comparing to COSMIC30 or the sum of two COSMIC30 sigs.
#'
#' @param sigs matrix with 96-channel substitution signatures as columns
#' @param threshold cosine similarity threshold
#' @return the list of signatures identified
#' @export
findClosestCOSMIC30andCombinations <- function(sigs,threshold){
  #load COSMIC30
  # cosmic30 <- read.table("../data/COSMIC_signatures.txt", sep="\t", header=T, as.is=T, check.names = FALSE)
  ncols30 <- length(cosmic30)
  colnames(cosmic30) <- 1:ncols30
  for(i in 1:(ncols30-1)){
    for(j in (i+1):ncols30){
      #message(paste0(i,"+",j))
      cosmic30[,paste0(i,"+",j)] <- (cosmic30[,i]+cosmic30[,j])/sum(cosmic30[,i]+cosmic30[,j])
    }
  }
  #compute cos sim matrix
  cos_sim_df <- data.frame()
  for (s in colnames(sigs)){
    for(a in colnames(cosmic30)){
      cos_sim_df[s,a] <- cos_sim(sigs[,s],cosmic30[,a])
    }
  }
  max.sim <- apply(cos_sim_df,1,max)
  closestCosmic <- colnames(cosmic30)[apply(cos_sim_df,1,which.max)]
  notFound <- max.sim < threshold
  closestCosmic[notFound] <- paste0("N",1:sum(notFound))
  closestCosmic[!notFound] <- paste0("C",closestCosmic[!notFound])
  return(closestCosmic)
}


#' Find closest COSMIC30 signatures
#'
#' Compares a set of signatures to the COSMIC30 signatures and
#' returns the list of signatures identified and the corresponding similarity. For example,
#' list(cosmic = c("C1","C3","C13"),cos_sim = c(0.94,0.85,0.7))
#' means that Cosmic (C) signatures 1, 3 and 13 were found, while
#' the corrsponding similarities to those signatures are 0.94, 0.85 and 0.7
#'
#' @param sigs matrix with 96-channel substitution signatures as columns
#' @return the list of signatures identified and corresponding similarities
#' @export
findClosestCOSMIC30_withSimilarity <- function(sigs){
  #load COSMIC30
  # cosmic30 <- read.table("../data/COSMIC_signatures.txt", sep="\t", header=T, as.is=T, check.names = FALSE)
  #compute cos sim matrix
  cos_sim_df <- data.frame()
  for (s in colnames(sigs)){
    for(a in colnames(cosmic30)){
      cos_sim_df[s,a] <- cos_sim(sigs[,s],cosmic30[,a])
    }
  }
  max.sim <- apply(cos_sim_df,1,max)
  closestCosmic <- apply(cos_sim_df,1,which.max)
  closestCosmic <- paste0("C",closestCosmic)
  res <- list()
  res[["cosmic"]] <- closestCosmic
  res[["cos.sim"]] <- max.sim
  return(res)
}

#' Find closest COSMIC30 signatures or combination of COSMIC30
#'
#' Compares a set of signatures to the COSMIC30 signatures or the simple sum of all combinations of two COSMIC signatures and
#' returns the list of signatures identified and the corresponding similarity. For example,
#' list(cosmic = c("C1","C3","C2+C13"),cos_sim = c(0.94,0.85,0.7))
#' means that Cosmic (C) signatures 1, 3 and 2+13 were found, while
#' the corrsponding similarities to those signatures are 0.94, 0.85 and 0.7
#'
#' @param sigs matrix with 96-channel substitution signatures as columns
#' @return the list of signatures identified and corresponding similarities
#' @export
findClosestCOSMIC30andCombinations_withSimilarity <- function(sigs){
  #load COSMIC30
  # cosmic30 <- read.table("../data/COSMIC_signatures.txt", sep="\t", header=T, as.is=T, check.names = FALSE)
  ncols30 <- length(cosmic30)
  colnames(cosmic30) <- 1:ncols30
  for(i in 1:(ncols30-1)){
    for(j in (i+1):ncols30){
      #message(paste0(i,"+",j))
      cosmic30[,paste0(i,"+",j)] <- (cosmic30[,i]+cosmic30[,j])/sum(cosmic30[,i]+cosmic30[,j])
    }
  }
  #compute cos sim matrix
  cos_sim_df <- data.frame()
  for (s in colnames(sigs)){
    for(a in colnames(cosmic30)){
      cos_sim_df[s,a] <- cos_sim(sigs[,s],cosmic30[,a])
    }
  }
  max.sim <- apply(cos_sim_df,1,max)
  closestCosmic <- colnames(cosmic30)[apply(cos_sim_df,1,which.max)]
  closestCosmic <- paste0("C",closestCosmic)
  res <- list()
  res[["cosmic"]] <- closestCosmic
  res[["cos.sim"]] <- max.sim
  return(res)
}

#returns the list of signatures identified. For example,
#c("R1","R3","N1","R5","N2")
#means that Rearrangement (R) signatures 1, 3 and 5 were found, while
#signatures N1 and N2 are unknown signatures (N for not found), based on the fact
#that no similarity >=threshold was found with the Rearr Sigs from Breast 560 study
findClosestRearrSigsBreast560 <- function(sigs,threshold){
  #load RS.Breast560
  # RS.Breast560 <- read.table("../data/rearrangement.signatures.txt", sep="\t", header=T, as.is=T, check.names = FALSE)
  #compute cos sim matrix
  cos_sim_df <- data.frame()
  for (s in colnames(sigs)){
    for(a in colnames(RS.Breast560)){
      cos_sim_df[s,a] <- cos_sim(sigs[,s],RS.Breast560[,a])
    }
  }
  max.sim <- apply(cos_sim_df,1,max)
  closestRS.Breast560 <- apply(cos_sim_df,1,which.max)
  notFound <- max.sim < threshold
  closestRS.Breast560[notFound] <- paste0("N",1:sum(notFound))
  closestRS.Breast560[!notFound] <- paste0("R",closestRS.Breast560[!notFound])
  return(closestRS.Breast560)
}

#returns the list of signatures identified and the corresponding similarity. For example,
#res$RS.Breast560 = c("R1","R3","R5")
#res$cos.sim = c(0.94,0.85,0.7)
#means that Rearrangement (R) signatures 1, 3 and 5 were found, while
#the corrsponding similarities to those signatures are 0.94, 0.85 and 0.7
findClosestRearrSigsBreast560_withSimilarity <- function(sigs){
  #load RS.Breast560
  # RS.Breast560 <- read.table("../data/rearrangement.signatures.txt", sep="\t", header=T, as.is=T, check.names = FALSE)
  #compute cos sim matrix
  cos_sim_df <- data.frame()
  for (s in colnames(sigs)){
    for(a in colnames(RS.Breast560)){
      cos_sim_df[s,a] <- cos_sim(sigs[,s],RS.Breast560[,a])
    }
  }
  max.sim <- apply(cos_sim_df,1,max)
  closestRS.Breast560 <- apply(cos_sim_df,1,which.max)
  closestRS.Breast560 <- paste0("R",closestRS.Breast560)
  res <- list()
  res[["RS.Breast560"]] <- closestRS.Breast560
  res[["cos.sim"]] <- max.sim
  return(res)
}

#' KL-divergence
#'
#' Compute the Kullback-Leibler Divergence between two matrices. In order to compute the divergence, .Machine$double.eps is added to matrices zero entries.
#'
#' @param m1 original matrix
#' @param m2 matrix to be used to approximate m1
#' @return KL-Divergence
#' @export
KLD <- function(m1,m2){
  # print(sessionInfo())
  # print(m1)
  # print(m2)
  # m1 <- as.vector(as.matrix(cat))
  # m2 <- as.vector(as.matrix(m2))
  m1[m1==0] <- .Machine$double.eps
  m2[m2==0] <- .Machine$double.eps
  return(sum(m1*(log(m1)-log(m2)) - m1 + m2))
}

#samples/sigantures are ararnged by columns

#' computeCorrelationOfTwoSetsOfSigs
#'
#' Compute the cosine similarity between two sets of signatures, which results in a cosine similarity matrix.
#'
#' @param sig1 matrix of signatures, with signatures as columns
#' @param sig2 matrix of signatures, with signatures as columns
#' @return cosine similarity matrix
#' @export
computeCorrelationOfTwoSetsOfSigs <- function(sigs1,sigs2){
  cos_sim_df <- data.frame()
  for (s in colnames(sigs1)){
    for(a in colnames(sigs2)){
      cos_sim_df[s,a] <- cos_sim(sigs1[,s],sigs2[,a])
    }
  }
  return(cos_sim_df)
}

removeSimilarCatalogueSamples <- function(cat,cosSimThreshold = 0.99){
  catCosSimCorr <- computeCorrelation(cat) - diag(ncol(cat))
  newcat <- data.frame(cat[,1],row.names = row.names(cat))
  colnames(newcat) = colnames(cat)[1]
  #as you are building the new catalogue, add samples only if they have cos sim
  #less than cosSimThreshold w.r.t. samples in the new catalgue
  for (j in 2:ncol(cat)){
    cossimres <- catCosSimCorr[j,colnames(newcat)]
    if(!any(cossimres>cosSimThreshold)){
      newcat <- cbind(newcat,cat[,j])
      colnames(newcat)[ncol(newcat)] <- colnames(cat)[j]
    }
  }
  return(newcat)
}

normaliseSamples <- function(cat){
  return(cat/matrix(rep(apply(cat,2,sum),nrow(cat)),nrow = nrow(cat),byrow = TRUE))
}

#' RMSE
#'
#' Function to compute the root mean squared error between two matrices.
#'
#' @param m1 first matrix to compare
#' @param m2 second matrix to compare
#' @return root mean squared error
#' @export
RMSE <- function(m1,m2){
  sqrt(sum((m1-m2)^2)/(ncol(m1)*nrow(m1)))
}

#' writeTable
#'
#' Utility function for simple write table with the following parameters: (sep = "\\t",quote = FALSE,row.names = TRUE,col.names = TRUE).
#'
#' @param t R table or matrix
#' @param file name of the output plain text file
#' @export
writeTable <- function(t,file){
  write.table(t,file = file,sep = "\t",quote = FALSE,row.names = TRUE,col.names = TRUE)
}

#' readTable
#'
#' Utility function for simple read table with the following parameters: (sep = "\\t",check.names = FALSE,header = TRUE,stringsAsFactors = FALSE).
#'
#' @param file name of the plain text file to read
#' @export
readTable <- function(file){
  read.table(file = file,sep = "\t",check.names = FALSE,header = TRUE,stringsAsFactors = FALSE)
}

#function to add shadowtext.
#Copied from https://github.com/cran/TeachingDemos/blob/master/R/shadowtext.R
#author: Greg Snow <538280@gmail.com>
shadowtext <- function(x, y=NULL, labels, col='white', bg='black',
                       theta= seq(pi/4, 2*pi, length.out=8), r=0.1, ... ) {

  xy <- xy.coords(x,y)
  xo <- r*strwidth('A')
  yo <- r*strheight('A')

  for (i in theta) {
    text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
  }
  text(xy$x, xy$y, labels, col=col, ... )
}

#' getOrganSignatures
#'
#' This function returns the organ-specific signatures for a given organ and mutation type as defined in Degasperi et al. 2020 Nat Cancer paper.
#'
#' @param typemut either subs, DNV or rearr
#' @param organ one of the following: "Biliary", "Bladder", "Bone_SoftTissue", "Breast", "Cervix", "CNS", "Colorectal", "Esophagus", "Head_neck", "Kidney", "Liver", "Lung", "Lymphoid", "Ovary", "Pancreas", "Prostate", "Skin", "Stomach", "Uterus"
#' @param version version "1" includes subs or rearr (ICGC cohort) organ-specific signatures from Degasperi et al. 2020, while version "2" includes the improved subs organ-specific signatures from ICGC as well as Hartwig and GEL, and the new DNV signatures. Set to "latest" to get the latest signature available for a given mutation type.
#' @param cohort for version 1 signatures only ICGC cohort is available, while for version 2 signatures ICGC, Hartwig and GEL cohort can be requested. Use "best" to get the most appropriate cohort for a given organ.
#' @return organ-specific signatures matrix
#' @references A. Degasperi, T. D. Amarante, J. Czarnecki, S. Shooter, X. Zou, D. Glodzik, S. Morganella, A. S. Nanda, C. Badja, G. Koh, S. E. Momen, I. Georgakopoulos-Soares, J. M. L. Dias, J. Young, Y. Memari, H. Davies, S. Nik-Zainal. A practical framework and online tool for mutational signature analyses show intertissue variation and driver dependencies, Nature Cancer, https://doi.org/10.1038/s43018-020-0027-5, 2020.
#' @export
getOrganSignatures <- function(organ,version="latest",cohort="best",typemut="subs",verbose = TRUE){
  sigs <- NULL
  if(typemut=="subs" & version=="1" & (cohort=="best" | cohort=="ICGC")){
    sigs <- all_organ_sigs_subs[,colnames(all_organ_sigs_subs)[grep(pattern = paste0("^",organ),colnames(all_organ_sigs_subs))]]
  }else if(typemut=="rearr" & (version=="1" | version=="latest") & (cohort=="best" | cohort=="ICGC")){
    sigs <- all_organ_sigs_rearr[,colnames(all_organ_sigs_rearr)[grep(pattern = paste0("^",organ),colnames(all_organ_sigs_rearr))]]
  }else if(typemut=="DNV" & (version=="2" | version=="latest")){
    if(cohort=="best") cohort <- "GEL"
    sigs <- organSignaturesDBSv1.01[,grepl(colnames(organSignaturesDBSv1.01),pattern = paste0(cohort,"-",organ)),drop=F]
  }else if(typemut=="subs" & (version=="2" | version=="latest")){
    if(cohort=="best") {
      cohort <- "GEL"
      if(organ=="Esophagus" | organ=="Head_neck") cohort <- "ICGC"
    }
    sigs <- organSignaturesSBSv2.03[,grepl(colnames(organSignaturesSBSv2.03),pattern = paste0(cohort,"-",organ)),drop=F]
  }
  if(ncol(sigs)==0 & verbose) message("[warning getOrganSignatures] Organ ",organ, " not available for mutation type ",typemut, ", version ",version, " and cohort ",cohort,".")
  return(sigs)
}

#' convertExposuresFromOrganToRefSigs
#'
#' This function converts the exposures matrix obatined from fitting organ-specific signatures into reference signatures exposures.
#' The function will detect the version of the signatures automatically.
#'
#' @param typemut either subs, DNV or rearr
#' @param expMatrix exposures matrix obatined from fitting organ-specific signatures
#' @return exposure matrix converted in reference signatures exposures
#' @references A. Degasperi, T. D. Amarante, J. Czarnecki, S. Shooter, X. Zou, D. Glodzik, S. Morganella, A. S. Nanda, C. Badja, G. Koh, S. E. Momen, I. Georgakopoulos-Soares, J. M. L. Dias, J. Young, Y. Memari, H. Davies, S. Nik-Zainal. A practical framework and online tool for mutational signature analyses show intertissue variation and driver dependencies, Nature Cancer, https://doi.org/10.1038/s43018-020-0027-5, 2020.
#' @export
convertExposuresFromOrganToRefSigs <- function(expMatrix,typemut="subs"){
  exposures <- NULL
  if(typemut=="subs"){
    if(all(rownames(expMatrix) %in% rownames(conversion_matrix_subs))){
      exposures <- t(conversion_matrix_subs[rownames(expMatrix),]) %*% as.matrix(expMatrix)
    }else if(all(rownames(expMatrix) %in% rownames(conversionMatrixSBSv2.03))){
      exposures <- t(conversionMatrixSBSv2.03[rownames(expMatrix),]) %*% as.matrix(expMatrix)
    }
  }else if(typemut=="rearr"){
    exposures <- t(conversion_matrix_rearr[rownames(expMatrix),]) %*% as.matrix(expMatrix)
  }else if(typemut=="DNV"){
    exposures <- t(conversionMatrixDBSv1.01[rownames(expMatrix),]) %*% as.matrix(expMatrix)
  }
  return(exposures)
}

