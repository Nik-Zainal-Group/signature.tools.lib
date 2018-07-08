prepare.rearr.catalogue_fromAnnotatedBedpe <- function(sv_bedpe) {
  
  catalogue.labels <- c('clustered_del_1-10Kb', 'clustered_del_10-100Kb', 'clustered_del_100Kb-1Mb', 'clustered_del_1Mb-10Mb', 'clustered_del_>10Mb', 'clustered_tds_1-10Kb', 'clustered_tds_10-100Kb', 'clustered_tds_100Kb-1Mb', 'clustered_tds_1Mb-10Mb', 'clustered_tds_>10Mb', 'clustered_inv_1-10Kb', 'clustered_inv_10-100Kb', 'clustered_inv_100Kb-1Mb', 'clustered_inv_1Mb-10Mb', 'clustered_inv_>10Mb', 'clustered_trans', 'non-clustered_del_1-10Kb', 'non-clustered_del_10-100Kb', 'non-clustered_del_100Kb-1Mb', 'non-clustered_del_1Mb-10Mb', 'non-clustered_del_>10Mb', 'non-clustered_tds_1-10Kb', 'non-clustered_tds_10-100Kb', 'non-clustered_tds_100Kb-1Mb', 'non-clustered_tds_1Mb-10Mb', 'non-clustered_tds_>10Mb', 'non-clustered_inv_1-10Kb', 'non-clustered_inv_10-100Kb', 'non-clustered_inv_100Kb-1Mb', 'non-clustered_inv_1Mb-10Mb', 'non-clustered_inv_>10Mb', 'non-clustered_trans')
  
  all_catalogues <- as.data.frame(matrix(nrow = length(catalogue.labels),ncol = 0))
  rownames(all_catalogues) <- catalogue.labels
  
  for (sample_name in unique(sv_bedpe$sample)){
    sample.rearrs <- sv_bedpe[sv_bedpe$sample==sample_name,]
    
    rearr_catalogue <- as.data.frame(matrix(0,nrow = length(catalogue.labels),ncol = 1))
  
    if (nrow(sample.rearrs)>0) {
  
      label1 <- rep('non-clustered', nrow(sample.rearrs))
      label1[ sample.rearrs$is.clustered] <- 'clustered'
  
      label2 <- rep('', nrow(sample.rearrs))
      label2[ sample.rearrs$svclass=='deletion'] <- '_del'
      label2[ sample.rearrs$svclass=='translocation'] <- '_trans'
      label2[ sample.rearrs$svclass=='inversion'] <- '_inv'
      label2[ sample.rearrs$svclass=='tandem-duplication'] <- '_tds'
  
      label3 <- rep('', nrow(sample.rearrs))
      sample.rearrs$bkdist <- abs(sample.rearrs$start2 - sample.rearrs$start1)
      label3[ sample.rearrs$svclass!='translocation' & sample.rearrs$bkdist<=1e4] <- '_1-10Kb'
      label3[ sample.rearrs$svclass!='translocation' & sample.rearrs$bkdist>1e4 & sample.rearrs$bkdist<=1e5 ] <- '_10-100Kb'
      label3[ sample.rearrs$svclass!='translocation' & sample.rearrs$bkdist>1e5 & sample.rearrs$bkdist<=1e6 ] <- '_100Kb-1Mb'
      label3[ sample.rearrs$svclass!='translocation' & sample.rearrs$bkdist>1e6 & sample.rearrs$bkdist<=1e7 ] <- '_1Mb-10Mb'
      label3[ sample.rearrs$svclass!='translocation' & sample.rearrs$bkdist>1e7 ] <- '_>10Mb'
  
      sample.rearrs$catalogue.label <- paste0(label1, label2, label3)
  
      sample.table <- as.data.frame(table( sample.rearrs$catalogue.label),drop=FALSE)
      rownames(sample.table ) <- sample.table$Var
  
      rearr_catalogue <-  sample.table [as.character(catalogue.labels), 'Freq',drop=FALSE ]
  
    }
  
    rearr.catalogue <- rearr_catalogue
    rownames(rearr.catalogue) <- catalogue.labels
    colnames(rearr.catalogue) <- sample_name
    rearr.catalogue[is.na(rearr.catalogue)] <- 0
  
    all_catalogues <- cbind(all_catalogues,rearr.catalogue)
  }
  
  all_catalogues
}
