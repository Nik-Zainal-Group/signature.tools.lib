context("testing computation of HRDLOH from ascat file")

test_that("test ascatToHRDLOH() using a input file", {
  
  ascat_test.df <- read.table("ascat_test.csv", sep=",", header=F, as.is=T, check.names = FALSE)
  colnames(ascat_test.df) <- c('seg_no', 'Chromosome', 'chromStart', 'chromEnd', 'total.copy.number.inNormal', 'minor.copy.number.inNormal', 'total.copy.number.inTumour', 'minor.copy.number.inTumour')
  
  res <- signature.tools.lib::ascatToHRDLOH(ascat_test.df,"testSample")
  
  expected_res <- 19
  names(expected_res) <- "testSample"
  
  expect_equal( res, expected_res )
  
})
