


#' Clustering of Catalogues or Signatures
#' 
#' Cluster sample catalogues or signatures using hierarchical clustering with average linkage. The distance metric used is 1 - cosine similarity, so the clustering is based on the shape/direction of the catalogue/signature vectors.
#' 
#' @param samplescatalogues catalogues or signatures, with colnames as sample or signature names and rownames as the channel names
#' @param nclusters number of clusters, can be a single positive integer or a list of integers
#' @param outdir if specified, plot results into the outdir folder
#' @return returns a list containing clustering results for each number of clusters requested
#' @keywords clustering
#' @export
#' @examples
#' resClust <- cataloguesClustering(catalogues,2:10,"results/")
cataloguesClustering <- function(samplescatalogues,
                                 nclusters, #nclusters can be both a value and a vector
                                 outdir=NULL){
  # clustering
  dM <- 1 - computeCorrelation_parallel(samplescatalogues)
  
  h <- hclust(as.dist(dM),method = "average")
  
  # now I need to run for the requested clusters, better clean things up
  selected_nclusters <- nclusters > 0 & nclusters <= ncol(dM)
  if(!all(selected_nclusters)){
    message("[warning cataloguesClustering] Following requested nclusters values will be ignored: ",
            paste(nclusters[!selected_nclusters],collapse = ", "),". Values of nclusters need to be between 1 and the number of columns of samplescatalogues.")
    nclusters <- nclusters[selected_nclusters]
  }
  
  if(length(nclusters)==0 | !is.numeric(nclusters)) {
    message("[error cataloguesClustering] No valid values for nclusters found. Exit.")
    return(NULL)
  }
  # order them just in case
  nclusters <- nclusters[order(nclusters)]
  
  # prepare data structures
  clusters_table <- data.frame(row.names = colnames(samplescatalogues),stringsAsFactors = F)
  sw_table <- data.frame(row.names = as.character(nclusters))
  clusters_stats <- list()
  mean_catalogues <- list()
  
  for (i in nclusters) {
    cl <- cutree(h,i)
    
    clusters_table[names(cl),as.character(i)] <- cl
    sil <- cluster::silhouette(cl,as.dist(dM))
    
    if(i==1){
      sw_table[as.character(i),"ASW"] <- NA
      sw_table[as.character(i),"n"] <- i
    }else if(!is.na(summary(sil)$avg.width)){
      sw_table[as.character(i),"ASW"] <- summary(sil)$avg.width
      sw_table[as.character(i),"n"] <- i
    }else{
      sw_table[as.character(i),"ASW"] <- NA
      sw_table[as.character(i),"n"] <- i
    }
    
    mean_catalogues[[i]] <- list()
    clusters_stats[[i]] <- data.frame(row.names = as.character(1:i))
    for(j in 1:i){
      tmpCatalogues <- samplescatalogues[,which(cl==j),drop=FALSE]
      clusters_stats[[i]][as.character(j),"nclsamples"] <- ncol(tmpCatalogues)
      clusters_stats[[i]][as.character(j),"cltotmuts"] <- sum(tmpCatalogues)
      clusters_stats[[i]][as.character(j),"clmeanmuts"] <- sum(tmpCatalogues)/ncol(tmpCatalogues)
      if(i>1) clusters_stats[[i]][as.character(j),"asw"] <- summary(sil)$clus.avg.widths[j]
      tmpCatalogues <- apply(tmpCatalogues,2,function(x) x/abs(sum(x)))
      mean_catalogues[[i]][[j]] <- matrix(apply(tmpCatalogues,1,mean),
                                          ncol = 1,nrow = nrow(tmpCatalogues),
                                          dimnames = list(rownames(tmpCatalogues),paste0("mean pattern cl",j)))
    }
    mean_catalogues[[i]] <- do.call(cbind,mean_catalogues[[i]])
  }
  
  # plot if requested
  if(!is.null(outdir)){
    pointsize <- 12
    dir.create(outdir,showWarnings = F,recursive = T)
    signature.tools.lib::writeTable(clusters_table,paste0(outdir,"/clusters_table.tsv"),row.names = T)
    signature.tools.lib::writeTable(sw_table,paste0(outdir,"/sw_table.tsv"),row.names = F)
    
    pdf(paste0(outdir,"/sw_table.pdf"),width = 6,height = 4.5,pointsize = pointsize)
    par(mar=c(5,5,3,1))
    plot(sw_table$n,sw_table$ASW,type="l",lwd=2,xlab="nclusters",ylab = "average silhouette width",las=1)
    grid()
    lines(sw_table$n,sw_table$ASW,type="l",lwd=2)
    dev.off()
    
    maxlengthlable <- max(strwidth(colnames(samplescatalogues),units = "inch",ps = par(ps=pointsize)))
    hcd <- as.dendrogram(h)
    cairo_pdf(paste0(outdir,"clusteringDendrogram.pdf"),height = 4+maxlengthlable,width = 0.18*length(h$order)+1,pointsize = pointsize)
    par(mai=c(0.2+maxlengthlable,0.8,0.4,0.2))
    plot(hcd,xlab = "",ylab = "1 - cosine similarity",sub = "",cex.lab = 1,cex.axis = 1)
    dev.off()
    
    for (i in nclusters) {
      idir <- paste0(outdir,"/nclusters_",i,"/")
      dir.create(idir,showWarnings = F,recursive = T)
      writeTable(clusters_stats[[i]],paste0(idir,"clusters_stats.tsv"),row.names = F)
      writeTable(mean_catalogues[[i]],paste0(idir,"meanOfClusters.tsv"),row.names = F)
      plotSignatures(signature_data_matrix = mean_catalogues[[i]],
                     output_file = paste0(idir,"meanOfClusters.pdf"),
                     plot_sum = FALSE,ncolumns = 3,
                     add_to_titles = paste0("n=",clusters_stats[[i]]$nclsamples," meanmuts=",sprintf("%.0f",clusters_stats[[i]]$clmeanmuts)))
      
    }
    
  }
  
  res <- list()
  res$mean_catalogues <- mean_catalogues
  res$clusters_table <- clusters_table
  res$sw_table <- sw_table
  res$clusters_stats <- clusters_stats
  res$nclusters <- nclusters
  res$samplescatalogues <- samplescatalogues
  return(res)
}



