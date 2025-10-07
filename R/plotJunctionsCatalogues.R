#' Plot Junctions Catalogues
#'
#' Function to plot one or more junctions catalogues.
#'
#' @param signature_data_matrix matrix of signatures, signatures as columns and channels as rows
#' @param output_file set output file, should end with ".jpg", "png" or ".pdf". If output_file==null, output will not be to a file, but will still run the plot functions. The option output_file==null can be used to add this plot to a larger output file.
#' @param plot_sum whether the sum of the channels should be plotted. If plotting signatures this should be FALSE, but if plotting sample catalogues, this can be set to TRUE to display the number of mutations in each sample.
#' @param overall_title set the overall title of the plot
#' @param mar set the option par(mar=mar)
#' @param howManyInOnePage how many signatures or catalogues should be plotted on one page. Multiple pages are plotted if more signatures/catalogues to plot have been requested
#' @param ncolumns how many columns should be used to arrange the signatures/catalogues to plot
#' @export
plotJunctionsCatalogues <-function(signature_data_matrix,
                                   output_file = NULL,
                                   plot_sum = TRUE,
                                   overall_title = "",
                                   add_to_titles = NULL,
                                   mar=NULL,
                                   howManyInOnePage=100,
                                   ncolumns=1,
                                   textscaling = 1){
  #This function plots a set of signatures in a single file, three signatures for each row.
  #signature_data_matrix is a data frame that contains the rearrangement signatures.
  #                      The columns are the signatures, while the rows are the 32 features
  # colnames(signature_data_matrix) <- sapply(colnames(signature_data_matrix),function(x) if (nchar(x)>30) paste0(substr(x,1,22),"...") else x)
  if(!is.null(output_file)) plottype <- substr(output_file,nchar(output_file)-2,nchar(output_file))
  # nt_col = rgb(228,26,28, maxColorValue = 255)
  # mh_col = rgb(77,175,74, maxColorValue =255)
  # other_col = rgb(152,78,163, maxColorValue =255)
  
  # nt_col = "#F3C300"
  # mh_col = "#F38400"
  # other_col = "#C2B280"
  
  nt_col = "#875692"
  mh_col = "#BE0032"
  other_col = "#C2B280"
  
  # nt_col = "#A1CAF1"
  # mh_col = "#0067A5"
  # other_col = "#C2B280"
  non_clust_col = rgb(240,240,240, maxColorValue =255)
  rearr.colours <- rep(c(rep(nt_col,3),other_col,rep(mh_col,3)),2)
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
      }else if(plottype=="png"){
        png(output_file,width = ncolumns*800,height = nplotrows*500,res = 220)
      }else if(plottype=="pdf"){
        pdf(output_file,width = ncolumns*8,height = nplotrows*5+0.5,pointsize = 26)
      }
      par(mfrow = c(nplotrows, ncolumns),omi=c(0,0,0.5,0),cex=0.7)
    }
    sizes <- c("1-3",
               "4-10",
               ">10")
    sizes_names <- rep(c(rev(sizes),"",sizes),2)
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
                     cex.axis = textscaling,
                     col=rearr.colours,
                     beside = FALSE,
                     cex.names = 0.8*textscaling,
                     cex.main = 0.9*textscaling,
                     las = 2,
                     border = 0,
                     space = 0.1)

      #save old plot coordinates
      op <- par("usr")
      #set new coordinates
      par(usr = c(0, 1, 0, 1))
      #add graphics
      par(xpd=TRUE)
      start1 <- 0.035
      xsep = 0.202
      start1_text <- 0.11
      tr_size <- 0.059

      stop <- start1
      for(i in 1:2){
        start <- stop
        stop <- start + xsep
        rect(start, -0.14, stop, -0.02,col = nt_col,lwd = 0,border = NA)
        text(x = start+0.5*xsep,y = -0.08,"nti",col = "white",cex = textscaling)
        start <- stop
        stop <- start + tr_size
        rect(start, -0.14, stop, -0.02,col = other_col,lwd = 0,border = NA)
        text(x = start+0.5*tr_size,y = -0.08,"o",col = "white",cex = textscaling)
        start <- stop
        stop <- start + xsep
        rect(start, -0.14, stop, -0.02,col = mh_col,lwd = 0,border = NA)
        text(x = start+0.5*xsep,y = -0.08,"mhd",col = "white",cex = textscaling)
      }
      xsep2 <- 2*xsep+tr_size
      rect(start1, -0.26, start1+xsep2, -0.14,col = "black",lwd = 0,border = NA)
      text(x = start1+0.5*xsep2,y = -0.2,"clustered",col = "white",cex = textscaling)
      rect(start1+xsep2, -0.26, start1+2*xsep2, -0.14,col = non_clust_col,lwd = 0,border = NA)
      text(x = start1+1.5*xsep2,y = -0.2,"non-clustered",col = "black",cex = textscaling)
      
      # plot text
      tstep <- xsep/3
      tstart <- start1+tstep/2
      tpos <- seq(from = tstart,by = tstep,length.out = length(sizes_names))
      text(x = tpos,y = -0.3,adj = 1,labels = sizes_names,cex = 0.8*textscaling,srt = 90)
      
      #restore old coordinates
      par(usr = op)
    }
    title(main = overall_title,outer = TRUE,cex.main = 1.5*textscaling)
    if(!is.null(output_file)) dev.off()
  }
}
