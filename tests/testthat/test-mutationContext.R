context("testing find mutation context of SNVs")

test_that("test findContextSNV()", {
  
  testSNVdf <- read.table(file = "test_hrdetect_1/test_hrdetect_1.snv.simple.txt",
                          sep = "\t",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
  
  # expected_cat <- read.table(file="test.snv.tab",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  
  res_context <- findContextSNV(snv_table = testSNVdf,
                                context_length = 2,
                                mtype = "C>G",
                                genomev = "hg19")
  
  expect_equal(res_context$C,c(57,53,230,63,39))
  
})