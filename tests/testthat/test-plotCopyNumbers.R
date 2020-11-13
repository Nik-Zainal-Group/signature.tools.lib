test_that("test ascatToHRDLOH() and plotCopyNumbers", {
  
  ascat_test.df <- read.table("ascat_test.csv", sep=",", header=F, as.is=T, check.names = FALSE)
  colnames(ascat_test.df) <- c('seg_no', 'Chromosome', 'chromStart', 'chromEnd', 'total.copy.number.inNormal', 'minor.copy.number.inNormal', 'total.copy.number.inTumour', 'minor.copy.number.inTumour')
  
  res <- ascatToHRDLOH(ascat_test.df,"testSample",return.loc = TRUE)
  
  plotCopyNumbers(sv_df = ascat_test.df,
                  sample_name = "testSample",
                  filename = "test.plotCopyNumbers.pdf",
                  highlightRegions = res,highlightText = paste0("HRDLOH\n",nrow(res)))
  
  file.remove("test.plotCopyNumbers.pdf")
  
  #if no error happen this code can be reached
  expect_true(TRUE)
  
})