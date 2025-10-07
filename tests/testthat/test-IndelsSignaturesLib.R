context("testing Indels catalogues")

test_that("test tabToIndelsCatalogue89() using an input file, hg19", {

  indels <- read.table("test.tabindels.tsv",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  indels$Sample <- "testSample"
  
  res <- tabToIndelsCatalogue89(indels = indels,genome.v = "hg19")
  
  expected_catalogues <- read.table("test.tabindels_catalogue89.tsv",sep = "\t",header = TRUE,stringsAsFactors = FALSE)

  
  expect_equal(nrow(res$indels_unfiltered),50)
  expect_equal(res$catalogues, expected_catalogues)
  
})
