
#' Plot Indels Signatures or Catalogues
#'
#' Function to plot one or more indels signatures or catalogues.
#'
#' @param signature_data_matrix matrix of signatures or catalogues, signatures as columns and channels as rows
#' @param output_file set output file, should end with ".jpg" or ".pdf". If output_file==null, output will not be to a file, but will still run the plot functions. The option output_file==null can be used to add this plot to a larger output file.
#' @param plot_sum whether the sum of the channels should be plotted. If plotting signatures this should be FALSE, but if plotting sample catalogues, this can be set to TRUE to display the number of mutations in each sample.
#' @param overall_title set the overall title of the plot
#' @param add_to_titles text to be added to the titles of each catalogue plot
#' @param mar set the option par(mar=mar)
#' @param howManyInOnePage how many signatures or catalogues should be plotted on one page. Multiple pages are plotted if more signatures/catalogues to plot have been requested
#' @param ncolumns how many columns should be used to arrange the signatures/catalogues to plot
#' @export
plotIDSignatures <- function(signature_data_matrix,
                               output_file = NULL,
                               plot_sum = TRUE,
                               overall_title = "",
                               add_to_titles = NULL,
                               mar=NULL,
                               howManyInOnePage=100,
                               ncolumns=1){
  
  # Check that the matrix contains the correct Indels channels
  IDchannels <- c( "A[Ins(C):R0]A"       ,   "A[Ins(C):R0]T"     ,     "Ins(C):R(0,3)"    ,      "Ins(C):R(4,6)"     ,     "Ins(C):R(7,9)"  ,       
                   "A[Ins(T):R(0,4)]A"   ,   "A[Ins(T):R(0,4)]C"  ,    "A[Ins(T):R(0,4)]G"  ,    "C[Ins(T):R(0,4)]A"  ,    "C[Ins(T):R(0,4)]C" ,    
                   "C[Ins(T):R(0,4)]G"  ,    "G[Ins(T):R(0,4)]A"  ,    "G[Ins(T):R(0,4)]C"  ,    "G[Ins(T):R(0,4)]G"  ,    "A[Ins(T):R(5,7)]A" ,    
                   "A[Ins(T):R(5,7)]C"  ,    "A[Ins(T):R(5,7)]G"  ,    "C[Ins(T):R(5,7)]A"  ,    "C[Ins(T):R(5,7)]C"  ,    "C[Ins(T):R(5,7)]G" ,    
                   "G[Ins(T):R(5,7)]A"  ,    "G[Ins(T):R(5,7)]C"  ,    "G[Ins(T):R(5,7)]G"   ,   "A[Ins(T):R(8,9)]A"   ,   "A[Ins(T):R(8,9)]C" ,    
                   "A[Ins(T):R(8,9)]G"  ,    "C[Ins(T):R(8,9)]A"  ,    "C[Ins(T):R(8,9)]C"  ,    "C[Ins(T):R(8,9)]G"  ,    "G[Ins(T):R(8,9)]A" ,    
                   "G[Ins(T):R(8,9)]C"  ,    "G[Ins(T):R(8,9)]G"  ,    "Ins(2,4):R0"       ,     "Ins(5,):R0"      ,       "Ins(2,4):R1"     ,      
                   "Ins(5,):R1"         ,    "Ins(2,):R(2,4)"   ,      "Ins(2,):R(5,9)"    ,     "[Del(C):R1]A"     ,      "[Del(C):R1]T"    ,      
                   "[Del(C):R2]A"      ,     "[Del(C):R2]T"       ,    "[Del(C):R3]A"      ,     "[Del(C):R3]T"      ,     "[Del(C):R(4,5)]A" ,     
                   "[Del(C):R(4,5)]T"   ,    "[Del(C):R(1,5)]G"   ,    "Del(C):R(6,9)"     ,     "A[Del(T):R(1,4)]A"  ,    "A[Del(T):R(1,4)]C" ,    
                   "A[Del(T):R(1,4)]G"  ,    "C[Del(T):R(1,4)]A"  ,    "C[Del(T):R(1,4)]C"   ,   "C[Del(T):R(1,4)]G"  ,    "G[Del(T):R(1,4)]A"  ,   
                   "G[Del(T):R(1,4)]C"  ,    "G[Del(T):R(1,4)]G"  ,    "A[Del(T):R(5,7)]A"  ,    "A[Del(T):R(5,7)]C"  ,    "A[Del(T):R(5,7)]G"  ,   
                   "C[Del(T):R(5,7)]A"  ,    "C[Del(T):R(5,7)]C"  ,    "C[Del(T):R(5,7)]G"  ,    "G[Del(T):R(5,7)]A"  ,    "G[Del(T):R(5,7)]C"  ,   
                   "G[Del(T):R(5,7)]G" ,     "A[Del(T):R(8,9)]A"  ,    "A[Del(T):R(8,9)]C"  ,    "A[Del(T):R(8,9)]G"   ,   "C[Del(T):R(8,9)]A"  ,   
                   "C[Del(T):R(8,9)]C"  ,    "C[Del(T):R(8,9)]G"  ,    "G[Del(T):R(8,9)]A"  ,    "G[Del(T):R(8,9)]C"   ,   "G[Del(T):R(8,9)]G"  ,   
                   "Del(2,4):R1"       ,     "Del(5,):R1"        ,     "Del(2,8):U(1,2):R(2,4)", "Del(2,):U(1,2):R(5,9)",  "Del(3,):U(3,):R2"  ,    
                   "Del(3,):U(3,):R(3,9)" ,  "Del(2,5):M1"      ,      "Del(3,5):M2"       ,     "Del(4,5):M(3,4)"    ,    "Del(6,):M1"     ,       
                   "Del(6,):M2"         ,    "Del(6,):M3"      ,       "Del(6,):M(4,)"     ,     "Complex"  )
  if(nrow(signature_data_matrix)==length(IDchannels)){
    if(all(rownames(rownames(signature_data_matrix)) %in% IDchannels)) {
      # make sure the order is correct for plotting
      signature_data_matrix <- signature_data_matrix[IDchannels,,drop=F]
    }else{
      message("[error] wrong channels in plotIDSignatures")
      return(NULL)
    }
  }else{
    message("[error] wrong number of channels in plotIDSignatures")
    return(NULL)
  }
  
  plotcolours <- c(rgb(102,153,204,maxColorValue = 255),
                   rgb(238,204,101,maxColorValue = 255),
                   rgb(238,153,170,maxColorValue = 255),
                   rgb(0,68,136,maxColorValue = 255),
                   rgb(153,120,0,maxColorValue = 255),
                   rgb(238,51,119,maxColorValue = 255),
                   rgb(118,43,131,maxColorValue = 255),
                   rgb(0,0,0,maxColorValue = 255))
  if(!is.null(output_file)) plottype <- substr(output_file,nchar(output_file)-2,nchar(output_file))
  nchannelsForEachType <- c(5,27,6,10,27,6,7,1)
  channel.colours <- c(rep(plotcolours[1],nchannelsForEachType[1]),
                       rep(plotcolours[2],nchannelsForEachType[2]),
                       rep(plotcolours[3],nchannelsForEachType[3]),
                       rep(plotcolours[4],nchannelsForEachType[4]),
                       rep(plotcolours[5],nchannelsForEachType[5]),
                       rep(plotcolours[6],nchannelsForEachType[6]),
                       rep(plotcolours[7],nchannelsForEachType[7]),
                       rep(plotcolours[8],nchannelsForEachType[8]))
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
        cairo_pdf(output_file,width = ncolumns*8,height = nplotrows*3+0.5,pointsize = 26)
      }
      par(mfrow = c(nplotrows, ncolumns),omi=c(0,0,0.5,0),cex=0.5)
    }
    for (pos in 1:ncol(tmpmatrix)){
      if(is.null(mar)){
        par(mar=c(2.8,3.5,4,1))
      }else{
        par(mar=mar)
      }
      titlei <- colnames(tmpmatrix)[pos]
      if (!is.null(add_to_titles)) titlei <- paste0(titlei," ",tmpadd[pos])
      if (plot_sum) titlei <- paste0(titlei,"\n(",round(sum(tmpmatrix[,pos]))," Indels)")
      muttypes <- c("1bp C","1bp T","\u22652bp","1bp C","1bp T","\u22652bp","Mh","X")
      muttypesTextCol <- c(rep("black",3),rep("white",5))
      xlabels <- rep("",nrow(tmpmatrix))
      barplot(tmpmatrix[,pos],
              # main = titlei,
              # cex.main = 0.9,
              names.arg = xlabels,
              col=channel.colours,
              beside = TRUE,
              las=2,
              xaxs='i',
              cex.names = 1,border = NA,space = 0.2)
      title(main = titlei,cex.main = 0.9,line = 1.5)
      par(xpd=TRUE)
      par(usr = c(0, 1, 0, 1))
      recttop <- -0.015
      rectbottom <- -0.14
      # start1 <- 0.035
      # gap <- 0.155
      start1 <- 0
      gap <- 1/nrow(tmpmatrix) + 0.000012
      for(j in 1:length(muttypes)){
        rect(start1, rectbottom, start1+gap*nchannelsForEachType[j], recttop,col = plotcolours[j],border = NA)
        start1 <- start1+gap*nchannelsForEachType[j]
      }
      start1 <- 0
      for(j in 1:length(muttypes)){
        textpos <- start1+(gap*nchannelsForEachType[j])/2
        text(x = textpos,y = -0.08,labels = muttypes[j],col = muttypesTextCol[j],font = 2,cex = 0.9)
        start1 <- start1+gap*nchannelsForEachType[j]
      }
      toplabels <- c("Insertions","Deletions","X")
      toplablesTextCol <- c("black","white","white")
      toplablesBackCol <- c(rgb(221,221,221,maxColorValue = 255),
                            rgb(136,136,136,maxColorValue = 255),
                            rgb(0,0,0,maxColorValue = 255))
      toplablesNchannels <- c(38,50,1)
      start1 <- 0
      rectbottom <- 1.015
      recttop <- 1.14
      for(j in 1:length(toplabels)){
        rect(start1, rectbottom, start1+gap*toplablesNchannels[j], recttop,col = toplablesBackCol[j],border = NA)
        textpos <- start1+(gap*toplablesNchannels[j])/2
        text(x = textpos,y = 1.08,labels = toplabels[j],col = toplablesTextCol[j],font = 2,cex = 0.9)
        start1 <- start1+gap*toplablesNchannels[j]
      }
      # textposx <- 0.04+seq(8,88,16)/104
      # text(x = textposx[1:3],y = -0.09,labels = muttypes[1:3],col = "white",font = 2)
      # text(x = textposx[4:6],y = -0.09,labels = muttypes[4:6],col = "black",font = 2)
      #shadowtext(x = 0.04+seq(8,88,16)/104,y = rep(-0.09,6),labels = muttypes,col = "white",bg = "black",r=0.2)
      text(x=0.5,y=-0.25,labels = "Indel Types",font = 2,cex = 0.9)
      par(xpd=FALSE)
    }
    title(main = overall_title,outer = TRUE,cex.main = 1.5)
    if(!is.null(output_file)) dev.off()
  }
}
