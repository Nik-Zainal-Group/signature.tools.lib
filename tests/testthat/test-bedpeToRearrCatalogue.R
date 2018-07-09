context("testing function bedpeToRearrCatalogue()")

test_that("test bedpeToRearrCatalogue() using a simple input file", {
  
  sv_bedpe <- read.table("test.bedpe",sep = "\t",header = TRUE,
                         stringsAsFactors = FALSE,check.names = FALSE)
  expected_res <- read.table("test.cat",sep = "\t",header = TRUE,
                             stringsAsFactors = FALSE,check.names = FALSE)
  res <- signature.tools.lib::bedpeToRearrCatalogue(sv_bedpe)
  
  expect_equal( res, expected_res )
  
})
