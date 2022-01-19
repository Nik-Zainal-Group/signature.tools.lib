context("testing conversion from dataframe of SNV to SNV catalogue")

test_that("test tabToSNVcatalogue() to convert a tab file with SNV into an SNV catalogue.", {
  
  testSNVdf <- read.table(file = "test_hrdetect_1/test_hrdetect_1.snv.simple.txt",
                          sep = "\t",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
  
  expected_cat <- read.table(file="test.snv.tab",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  
  res_cat <- tabToSNVcatalogue(testSNVdf,genome.v = "hg19")
  
  expect_equal(res_cat$catalogue,expected_cat)
  
})