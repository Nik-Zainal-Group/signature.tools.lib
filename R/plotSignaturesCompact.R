#' Plot Substitution Signatures or Catalogues, Compact plot
#'
#' Function to plot one or more substitution signatures or catalogues. Minimal compact plot without labels and axis.
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
plotSubsSignaturesCompact <- function(signature_data_matrix,
                               output_file = NULL,
                               plot_sum = FALSE,
                               overall_title = "",
                               add_to_titles = NULL,
                               mar=NULL,
                               cex.title = 1,
                               cex.axistext = 1,
                               howManyInOnePage=100,
                               ncolumns=1){

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
        jpeg(output_file,width = ncolumns*600,height = nplotrows*280,res = 220)
      }else if(plottype=="pdf"){
        pdf(output_file,width = ncolumns*6,height = nplotrows*2.8,pointsize = 26)
      }
      par(mfrow = c(nplotrows, ncolumns),oma=c(0,0,0,0),cex=0.7)
    }
    for (pos in 1:ncol(tmpmatrix)){
      if(is.null(mar)){
        par(mar=c(1,0.5,1.5,0.5))
      }else{
        par(mar=mar)
      }
      title <- colnames(tmpmatrix)[pos]
      if (plot_sum) title <- paste0(title," (",round(sum(tmpmatrix[,pos]))," SNVs)")
      if (!is.null(add_to_titles)) title <- paste0(title,"\n",tmpadd[pos])
      muttypes <- c("C>A","C>G","C>T","T>A","T>C","T>G")
      xlabels <- rep("",96)
      barplot(tmpmatrix[,pos],
              main = title,
              names.arg = xlabels,
              col=rearr.colours,
              beside = TRUE,
              las=2,
              cex.names = 1,border = NA,space = 0.2,
              cex.axis = cex.axistext,cex.main = cex.title,yaxt='n')
      par(xpd=TRUE)
      par(usr = c(0, 1, 0, 1))
      recttop <- -0.02
      rectbottom <- -0.08
      start1 <- 0.035
      gap <- 0.155
      rect(start1, rectbottom, start1+gap, recttop,col = plotcolours[1],border = NA)
      rect(start1+gap, rectbottom, start1+2*gap, recttop,col = plotcolours[2],border = NA)
      rect(start1+2*gap, rectbottom, start1+3*gap, recttop,col = plotcolours[3],border = NA)
      rect(start1+3*gap, rectbottom, start1+4*gap, recttop,col = plotcolours[4],border = NA)
      rect(start1+4*gap, rectbottom, start1+5*gap, recttop,col = plotcolours[5],border = NA)
      rect(start1+5*gap, rectbottom, start1+6*gap, recttop,col = plotcolours[6],border = NA)
      textposx <- 0.04+seq(8,88,16)/104
      par(xpd=FALSE)
    }
    title(main = overall_title,outer = TRUE,cex.main = 1.5)
    if(!is.null(output_file)) dev.off()
  }
}

#' Plot Dinucleotide Variant Signatures or Catalogues, Compact plot
#'
#' Function to plot one or more DNV signatures or catalogues.  Minimal compact plot without labels and axis.
#'
#' @param signature_data_matrix matrix of signatures or catalogues, signatures as columns and channels as rows. You can use either Zou or Alexandrov Style and the function will recognise the rownames and plot accordingly
#' @param output_file set output file, should end with ".jpg" of ".pdf". If output_file==null, output will not be to a file, but will still run the plot functions. The option output_file==null can be used to add this plot to a larger output file.
#' @param plot_sum whether the sum of the channels should be plotted. If plotting signatures this should be FALSE, but if plotting sample catalogues, this can be set to TRUE to display the number of mutations in each sample.
#' @param overall_title set the overall title of the plot
#' @param add_to_titles text to be added to the titles of each catalogue plot
#' @param mar set the margin of the plot
#' @param howManyInOnePage how many signatures or catalogues should be plotted on one page. Multiple pages are plotted if more signatures/catalogues to plot have been requested
#' @param ncolumns how many columns should be used to arrange the signatures/catalogues to plot
#' @export
plotDNVSignaturesCompact <- function(signature_data_matrix,
                              output_file = NULL,
                              plot_sum = FALSE,
                              overall_title = "",
                              add_to_titles = NULL,
                              mar=NULL,
                              cex.title = 1,
                              cex.axistext = 1,
                              howManyInOnePage=100,
                              ncolumns=1){
  if(!is.null(output_file)) plottype <- substr(output_file,nchar(output_file)-2,nchar(output_file))

  mutationTypesZou <- c("AA>CC","AA>CG","AA>CT","AA>GC","AA>GG","AA>GT","AA>TC","AA>TG","AA>TT",
                        "AC>CA","AC>CG","AC>CT","AC>GA","AC>GG","AC>GT","AC>TA","AC>TG","AC>TT",
                        "AG>CA","AG>CC","AG>CT","AG>GA","AG>GC","AG>GT","AG>TA","AG>TC","AG>TT",
                        "AT>CA","AT>CC","AT>CG","AT>GA","AT>GC","AT>TA",
                        "CA>AC","CA>AG","CA>AT","CA>GC","CA>GG","CA>GT","CA>TC","CA>TG","CA>TT",
                        "CC>AA","CC>AG","CC>AT","CC>GA","CC>GG","CC>GT","CC>TA","CC>TG","CC>TT",
                        "CG>AA","CG>AC","CG>AT","CG>GA","CG>GC","CG>TA",
                        "GA>AC","GA>AG","GA>AT","GA>CC","GA>CG","GA>CT","GA>TC","GA>TG","GA>TT",
                        "GC>AA","GC>AG","GC>AT","GC>CA","GC>CG","GC>TA",
                        "TA>AC","TA>AG","TA>AT","TA>CC","TA>CG","TA>GC")
  mutationTypesAlexandrov <-    c("AC>CA","AC>CG","AC>CT","AC>GA","AC>GG","AC>GT","AC>TA","AC>TG","AC>TT",
                                  "AT>CA","AT>CC","AT>CG","AT>GA","AT>GC","AT>TA",
                                  "CC>AA","CC>AG","CC>AT","CC>GA","CC>GG","CC>GT","CC>TA","CC>TG","CC>TT",
                                  "CG>AT","CG>GC","CG>GT","CG>TA","CG>TC","CG>TT",
                                  "CT>AA","CT>AC","CT>AG","CT>GA","CT>GC","CT>GG","CT>TA","CT>TC","CT>TG",
                                  "GC>AA","GC>AG","GC>AT","GC>CA","GC>CG","GC>TA",
                                  "TA>AT","TA>CG","TA>CT","TA>GC","TA>GG","TA>GT",
                                  "TC>AA","TC>AG","TC>AT","TC>CA","TC>CG","TC>CT","TC>GA","TC>GG","TC>GT",
                                  "TG>AA","TG>AC","TG>AT","TG>CA","TG>CC","TG>CT","TG>GA","TG>GC","TG>GT",
                                  "TT>AA","TT>AC","TT>AG","TT>CA","TT>CC","TT>CG","TT>GA","TT>GC","TT>GG")

  if(all(rownames(signature_data_matrix)==mutationTypesZou)){
    catalogueType <- "Zou"
  }else if(all(rownames(signature_data_matrix)==mutationTypesAlexandrov)){
    catalogueType <- "Alexandrov"
  }else{
    message("Unknown DNV channels (rownames) style, please use NN>NN rownames in either Zou or Alexandrov style")
    exit(1)
  }

  if(catalogueType=="Zou"){
    muttypeslength <- c(9,9,9,6,9,9,6,9,6,6)
    mypalette <- c("#e6194b","#3cb44b","#ffe119","#0082c8","#f58231","#911eb4","#46f0f0","#f032e6","#d2f53c","#fabebe")
    muttypes <- c("AA>NN","AC>NN","AG>NN","AT>NN","CA>NN","CC>NN","CG>NN","GA>NN","GC>NN","TA>NN")
  }else if(catalogueType=="Alexandrov"){
    muttypeslength <- c(9,6,9,6,9,6,6,9,9,9)
    mypalette <- c(rgb(164,205,224,maxColorValue = 255),
                   rgb(31,119,182,maxColorValue = 255),
                   rgb(176,222,139,maxColorValue = 255),
                   rgb(50,160,43,maxColorValue = 255),
                   rgb(252,152,152,maxColorValue = 255),
                   rgb(227,33,29,maxColorValue = 255),
                   rgb(243,185,109,maxColorValue = 255),
                   rgb(250,125,0,maxColorValue = 255),
                   rgb(200,175,210,maxColorValue = 255),
                   rgb(105,60,155,maxColorValue = 255))
    muttypes <- c("AC>NN","AT>NN","CC>NN","CG>NN","CT>NN","GC>NN","TA>NN","TC>NN","TG>NN","TT>NN")
  }else{
    message("Unknown catalogue type.")
    exit(1)
  }

  rearr.colours <- c()
  for (i in 1:length(mypalette)) rearr.colours <- c(rearr.colours,rep(mypalette[i],muttypeslength[i]))
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
        jpeg(output_file,width = ncolumns*600,height = nplotrows*280,res = 220)
      }else if (plottype=="pdf"){
        pdf(output_file,width = ncolumns*6,height = nplotrows*2.8,pointsize = 26)
      }
      par(mfrow = c(nplotrows, ncolumns),oma=c(0,0,0,0),cex=0.7)
    }
    for (pos in 1:ncol(tmpmatrix)){
      if(is.null(mar)){
        par(mar=c(1,0.5,1.5,0.5))
      }else{
        par(mar=mar)
      }
      title <- colnames(tmpmatrix)[pos]
      if (plot_sum) title <- paste0(title," (",round(sum(tmpmatrix[,pos]))," DNVs)")
      if (!is.null(add_to_titles)) title <- paste0(title,"\n",tmpadd[pos])
      xlabels <- rep("",nrow(tmpmatrix))
      xlabels2 <- sapply(rownames(tmpmatrix),function(x){
        strsplit(x,split = ">")[[1]][2]
      })

      b <- barplot(tmpmatrix[,pos],
                   main = title,
                   names.arg = xlabels,
                   col=rearr.colours,
                   beside = TRUE,
                   las=2,
                   cex.names = 1,border = NA,space = 0.2,
                   cex.axis = cex.axistext,cex.main = cex.title,yaxt='n')
      par(xpd=TRUE)
      par(usr = c(0, 1, 0, 1))
      recttop <- -0.02
      rectbottom <- -0.08
      start1 <- 0.037
      endfinal <- 0.963
      gap <- (endfinal-start1)/nrow(tmpmatrix)
      xpos2 <- start1 - gap/2 + 1:length(xlabels2)*gap
      for (i in 1:length(mypalette)) {
        end1 <- start1+gap*muttypeslength[i]
        rect(start1, rectbottom, end1, recttop,col = mypalette[i],lwd = 0,border = NA)
        start1 <- end1
      }

      par(xpd=FALSE)
    }
    title(main = overall_title,outer = TRUE,cex.main = 1.5)
    if(!is.null(output_file)) dev.off()
  }
}
