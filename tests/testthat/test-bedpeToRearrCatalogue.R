context("testing building a rearrangement catalogue from a BEDPE file")

test_that("test bedpeToRearrCatalogue() using a simple input file", {
  
  sv_bedpe <- read.table("test.bedpe",sep = "\t",header = TRUE,
                         stringsAsFactors = FALSE,check.names = FALSE)
  expected_res <- read.table("test.cat",sep = "\t",header = TRUE,
                             stringsAsFactors = FALSE,check.names = FALSE)
  res <- signature.tools.lib::bedpeToRearrCatalogue(sv_bedpe)
  
  expect_equal( res, expected_res )
  
})
