context("testing conversion from dataframe of SNV to DNV catalogue")

test_that("test snvTabToDNVcatalogue() to convert a tab file with SNV into an DNV catalogue.", {
  
  testSNVdf <- read.table(file = "test_DNVs.snv.simple.txt",
                          sep = "\t",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
  
  expected_cat <- read.table(file="test.dnv.expected.tab",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  
  res_cat <- snvTabToDNVcatalogue(testSNVdf)
  
  expect_equal(res_cat$DNV_catalogue,expected_cat)
  
})