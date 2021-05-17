context("testing functions related to organ specific signatures")

test_that("test getOrganSignatures() latest SBS sigs", {
  
  os1 <- getOrganSignatures("Breast")
  expect_true(ncol(os1)>0)
  
})

test_that("test getOrganSignatures() RefSig v1 subs", {
  
  os1 <- getOrganSignatures("Breast",version = "1")
  expect_true(ncol(os1)>0)
  
})

test_that("test getOrganSignatures() latest DBS sigs", {
  
  os1 <- getOrganSignatures("Breast",typemut = "DNV")
  expect_true(ncol(os1)>0)
  
})

test_that("test getOrganSignatures() latest Rearr sigs", {
  
  os1 <- getOrganSignatures("Breast",typemut = "rearr")
  expect_true(ncol(os1)>0)
  
})

test_that("test convert organ-specific exposures to RefSig SBS exposures", {
  
  expo1 <- readTable("exposures/exposuresExampleOrganSpecificV2.tsv")
  expo1 <- convertExposuresFromOrganToRefSigs(expo1)
  expect_true(nrow(expo1)>0)
  
})
