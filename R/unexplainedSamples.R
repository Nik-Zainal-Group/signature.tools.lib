
#' Estimate samples not fully explained by the given signatures
#'
#' Given a catalogue of samples, and a set of signatures, this functions identifies
#' samples that have a  significantly higher error than the rest. It is important to
#' have enough samples that are fully explained by the singatures in the catalogue
#' matrix, so that the background error distribution is well defined and the samples
#' with significantly higher error can be identified more clearly.
#' Two criteria are used to identify the samples that are not  fully explained by the
#' signatures. The first criterion considers the relative amount of mutations and is governed
#' by the parameters pvalueMethod and pvalue_threshold. In practice, we calculate a p-value
#' for each sample error or residual (normalised by total mutations of the sample) to
#' determine if i is higher than the rest. The parameter pvalueMethod allows to select
#' different errors, though we advise to use normErrorSAD, which is the sum of absolute
#' deviations of every channel divided by the number of mutations in the sample.
#' This first criterion can be disabled by setting considerOnlyNmutsThreshold=TRUE.
#' The second criterion is a minimum number of mutations that should be in the error
#' or residual, which can be set with nmuts_threshold. The parameter nmutsMethod
#' determines the type of error used, here we suggest residualSSD, the sum of signed
#' differences of every channel of the residual.
#' In this function we refer to the error as the difference between the
#' original catalogue and the catalogue reconstructed from the exposures obtained
#' with a simple signature fit (optimise for KLD), while we refer to the residual
#' as the difference between the original catalogue and the catalogue reconstructed
#' from the exposures obtained with a contrained fit, where the difference is
#' constrained to be mostly positive, to better highlight the pattern that might be
#' present in the catalogue butt not captured by he given signatures.
#' Both criteria are necessary to consider a sample unexplained, unless considerOnlyNmutsThreshold=TRUE.
#'
#' @param outfileRoot if specified, generate a plot, otherwise no plot is generated
#' @param catalogues original catalogues, channels as rows and samples as columns
#' @param sigs mutational signatures used for fitting, channels as rows, signatures as columns
#' @param nmuts_threshold minimum number of mutations in the error/residual to consider the samples unexplained
#' @param pvalue_threshold threshold for statistical significance of the normalised error/residual
#' @param pvalueMethod method to be used for the relative criterion. Default is normErrorSAD. Alternatives are normResidualSAD or normResidualSSD
#' @param nmutsMethod method to be used for the absolute criterion. Default is residualSSD. Alternatives are errorSAD or residualSAD
#' @param considerOnlyNmutsThreshold use only the absolute criterion and disables the relative criterion
#' @keywords unexplained samples
#' @return table of samples with associated error metrics and samples with significant error and/orr residual highlighted
#' @export
#' @examples
#' resObj <- unexplainedSamples(catalogues=catalogues,
#'                              sigs=signatures)
unexplainedSamples <- function(outfileRoot=NULL,
                               catalogues,
                               sigs,
                               pvalue_threshold=0.03,
                               nmuts_threshold=300,
                               pvalueMethod="normErrorSAD", #normErrorSAD or normResidualSAD or normResidualSSD
                               nmutsMethod="residualSSD", #residualSSD or errorSAD or residualSAD
                               considerOnlyNmutsThreshold=FALSE){
  # initial fit with common signatures
  exposures <- SignatureFit(catalogues,signature_data_matrix = sigs,method = "KLD",showDeprecated = F,doRound = F)
  # initial fit with common signatures using Tau's method
  exposures_constr <- flexconstr_sigfit_multipleSamples(as.matrix(sigs),as.matrix(catalogues),allmut_tolratio = 0.003)
  # calculate the residual for all samples
  residuals <- catalogues - as.matrix(sigs) %*% as.matrix(exposures_constr)
  # calculate the error from the reconstucted catalogues
  error <- catalogues - as.matrix(sigs) %*% as.matrix(exposures)

  #look at various metrics
  # SAD = Sum of Absolute Deviations
  # SSD = Sum of Signed Deviations
  # errorSSD tends to zero so it is useless here
  errorSAD <- c()
  residualSAD <- c()
  residualSSD <- c()
  norm_errorSAD <- c()
  norm_residualSAD <- c()
  norm_residualSSD <- c()
  for (i in 1:ncol(catalogues)){
    errorSAD <- c(errorSAD,sum(abs(error[,i])))
    residualSAD <- c(residualSAD,sum(abs(residuals[,i])))
    residualSSD <- c(residualSSD,sum(residuals[,i]))
    sumcati <- sum(catalogues[,i])
    norm_errorSAD <- c(norm_errorSAD,errorSAD[i]/sumcati)
    norm_residualSAD <- c(norm_residualSAD,residualSAD[i]/sumcati)
    norm_residualSSD <- c(norm_residualSSD,residualSSD[i]/sumcati)
  }

  # use requested values to determine the unexplained samples
  if(pvalueMethod=="normErrorSAD"){
    normValues <- norm_errorSAD
    xlab <- "-log10(p-value) of error SAD/nmuts"
  }else if(pvalueMethod=="normResidualSAD"){
    normValues <- norm_residualSAD
    xlab <- "-log10(p-value) of residual SAD/nmuts"
  }else if(pvalueMethod=="normResidualSSD"){
    normValues <- norm_residualSSD
    xlab <- "-log10(p-value) of residual SSD/nmuts"
  }else {
    message("Invalid pvalueMethod ",pvalueMethod,", plese use one of these: normErrorSAD, normResidualSAD or normResidualSSD")
    return(NULL)
  }

  if(nmutsMethod=="residualSSD"){
    values <- residualSSD
    ylab <- "residual SSD"
  }else if(nmutsMethod=="errorSAD"){
    values <- errorSAD
    ylab <- "error SAD"
  }else if(nmutsMethod=="residualSAD"){
    values <- residualSAD
    ylab <- "residual SAD"
  }else {
    message("Invalid nmutsMethod ",nmutsMethod,", plese use one of these: residualSSD, errorSAD or residualSAD")
    return(NULL)
  }

  # get p-values for the norm values
  pval_normValues <- c()
  for (i in 1:length(normValues)) {
    m <- mean(normValues[-i])
    s <- sd(normValues[-i])
    pval_normValues <- c(pval_normValues, 1 - pnorm(normValues[i],
                                                    mean = m, sd = s, lower.tail = TRUE))
  }

  selection <- values >= nmuts_threshold
  if(!considerOnlyNmutsThreshold) selection <- selection & pval_normValues < pvalue_threshold
  which_significant <- which(selection)


  #prepare return obj
  res <- list()
  res$info_samples <- data.frame(index=1:ncol(catalogues),
                                 sample=colnames(catalogues),
                                 isUnexplainedSample=selection,
                                 nmutsMethod=rep(nmutsMethod,ncol(catalogues)),
                                 pvalueMethod=rep(pvalueMethod,ncol(catalogues)),
                                 pval_normValues=pval_normValues,
                                 errorSAD=errorSAD,
                                 residualSAD=residualSAD,
                                 residualSSD=residualSSD,
                                 norm_errorSAD=norm_errorSAD,
                                 norm_residualSAD=norm_residualSAD,
                                 norm_residualSSD=norm_residualSSD,
                                 stringsAsFactors = F)
  res$which_significant <- which_significant
  res$exposures <- exposures
  res$unexplSamples <- colnames(catalogues)[which_significant]
  res$unexplSamples_catalogues <- catalogues[,which_significant,drop=F]
  res$unexplSamples_residuals <- residuals[,which_significant,drop=F]
  res$unexplSamples_errors <- error[,which_significant,drop=F]
  res$numberOfUnexplSamples <- length(which_significant)
  res$all_residuals <- residuals
  res$all_errors <- error

  # also plot if requested
  if(!is.null(outfileRoot)){

    trimdir <- function(x){
      while(substr(x,nchar(x),nchar(x))!="/" & nchar(x)>0){
        x <- substr(x,1,nchar(x)-1)
      }
      return(x)
    }

    # create dir
    xdir <- trimdir(outfileRoot)
    if(nchar(xdir)>0) dir.create(xdir,recursive = T,showWarnings = F)

    if(considerOnlyNmutsThreshold){
      plottitle <- bquote(atop(paste("Unexplained Catalogues"),paste(.(nmutsMethod)>=.(nmuts_threshold))))
    }else{
      plottitle <- bquote(atop(paste("Unexplained Catalogues"),paste(.(nmutsMethod)>=.(nmuts_threshold),", ",.(pvalueMethod)," p-value"<.(pvalue_threshold),sep="")))
    }
    cairo_pdf(filename = paste0(outfileRoot,"_UnexplainedCataloguesSelectionPlot.pdf"),width = 7,height = 5,pointsize = 14)
    par(mar=c(4,6,3.5,2),mgp=c(2.5,1,0))
    plot(-log10(pval_normValues[values>0]),values[values>0],
         col=rgb(0.5,0.5,0.5,0.5),
         pch = 16,
         xlab = xlab,las=1,
         ylab = "",
         log="y",
         main = plottitle,cex.main=0.9)
    title(ylab = ylab,mgp=c(4,1,0))
    if(length(which_significant)>0) points(-log10(pval_normValues[which_significant]),values[which_significant],col="red",pch = 16)
    if(!considerOnlyNmutsThreshold) abline(v = -log10(pvalue_threshold),col = "red",lty = 2)
    abline(h = nmuts_threshold,col = "red",lty = 2)
    dev.off()

    plotError(outfileRoot = outfileRoot,
              catalogues = catalogues,
              sigs = sigs,
              samplesToHighlight = res$unexplSamples,
              exposures = exposures)

    if(length(which_significant)>0){
      plotSignatures(output_file = paste0(outfileRoot,"_UnexplainedCatalogues.pdf"),res$unexplSamples_catalogues,ncolumns = 3)
      plotSignatures(output_file = paste0(outfileRoot,"_UnexplainedResiduals.pdf"),res$unexplSamples_residuals,ncolumns = 3)
      plotSignatures(output_file = paste0(outfileRoot,"_UnexplainedError.pdf"),res$unexplSamples_errors,ncolumns = 3)
    }
    write.table(x = res$info_samples,
                file = paste0(outfileRoot,"_UnexplainedCataloguesSelection.tsv"),sep = "\t",
                quote = FALSE,col.names = T,row.names = F)
  }
  return(res)
}


plotError <- function(outfileRoot=NULL,
                      catalogues,
                      sigs,
                      exposures,
                      samplesToHighlight=NULL,
                      tag=""){
  error <- catalogues - as.matrix(sigs) %*% as.matrix(exposures)
  errorSAD <- c()
  norm_errorSAD <- c()
  for (i in 1:ncol(catalogues)){
    errorSAD <- c(errorSAD,sum(abs(error[,i])))
    sumcati <- sum(catalogues[,i])
    norm_errorSAD <- c(norm_errorSAD,errorSAD[i]/sumcati)
  }
  plottitle <- "Normalised Sum of Absolute Deviations (SAD)\nbetween original and reconstructed catalogue"
  xlab <- "samples"
  ylab <- "SAD/n. mutations"
  cairo_pdf(filename = paste0(outfileRoot,"_ErrorAcrossSamples",tag,".pdf"),width = 7,height = 5,pointsize = 14)
  par(mar=c(4,6,3.5,2),mgp=c(2.5,1,0))
  plot(1:length(norm_errorSAD),norm_errorSAD,
       col=rgb(0.5,0.5,0.5,0.5),
       pch = 16,
       xlab = xlab,las=1,
       ylab = "",
       log="y",
       main = plottitle,cex.main=0.9)
  title(ylab = ylab,mgp=c(4,1,0))
  if(!is.null(samplesToHighlight)){
    selectedPoints <- which(colnames(catalogues) %in% samplesToHighlight)
    points(selectedPoints,norm_errorSAD[selectedPoints],col="red",pch=16)
  }
  dev.off()

}

