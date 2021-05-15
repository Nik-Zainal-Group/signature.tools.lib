context("testing DNV and TNV functions")

test_that("test dnvTabToDNVcatalogue() to convert a tab file with DNV into an DNV catalogue.", {
  
  testDNVdf <- read.table(file = "DNVTNVexample/faff4626-615b-416a-b7a6-9d177dcc94a9_DNV.tsv",
                          sep = "\t",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
  
  expected_cat <- read.table(file="DNVTNVexample/faff4626-615b-416a-b7a6-9d177dcc94a9_DNVcatalogue.tsv",
                             sep = "\t",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
  
  res_cat <- dnvTabToDNVcatalogue(testDNVdf)
  
  expect_equal(res_cat$DNV_catalogue,expected_cat)
  
})

test_that("test dnvTabToTNVcatalogue() to convert a tab file with DNV into an TNV catalogue.", {
  
  testDNVdf <- read.table(file = "DNVTNVexample/faff4626-615b-416a-b7a6-9d177dcc94a9_DNV.tsv",
                          sep = "\t",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
  
  expected_cat <- read.table(file="DNVTNVexample/faff4626-615b-416a-b7a6-9d177dcc94a9_TNVcatalogue.tsv",
                             sep = "\t",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
  
  res_cat <- dnvTabToTNVcatalogue(testDNVdf)
  
  expect_equal(res_cat$TNV_catalogue,expected_cat)
  
})

test_that("test tnvTabToTNVcatalogue() to convert a tab file with TNV into an TNV catalogue.", {
  
  testTNVdf <- read.table(file = "DNVTNVexample/faff4626-615b-416a-b7a6-9d177dcc94a9_TNV.tsv",
                          sep = "\t",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
  
  expected_cat <- read.table(file="DNVTNVexample/faff4626-615b-416a-b7a6-9d177dcc94a9_TNVcatalogue.tsv",
                             sep = "\t",header = TRUE,stringsAsFactors = FALSE,check.names = FALSE)
  
  res_cat <- tnvTabToTNVcatalogue(testTNVdf)
  
  expect_equal(res_cat$TNV_catalogue,expected_cat)
  
})