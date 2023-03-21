context("testing Indels classification from data frame")

test_that("test tabToIndelsClassification() using an input file, check count and proportion of indels is correct, hg19", {

  indels <- read.table("test.tabindels.tsv",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  
  res <- tabToIndelsClassification(indels,"testSample","hg19")
  
  expect_count_proportion <- data.frame(sample = "testSample",
                                        del.mh = 22,
                                        del.rep = 15,
                                        del.none = 4,
                                        ins = 8,
                                        complex = 1,
                                        all.del = 41,
                                        all.indels = 50,
                                        del.mh.prop = 0.536585365853659,
                                        del.rep.prop = 0.365853658536585,
                                        del.none.prop = 0.0975609756097561,
                                        del.mh.count = 22,
                                        del.rep.count = 15,
                                        del.none.count = 4)

  
  
  expect_equal( res$count_proportion, expect_count_proportion )
  
})

test_that("test tabToIndelsClassification() using an input file, check count and proportion of indels is correct, hg38", {
  
  indels <- read.table("test.tabindels.tsv",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  indels$chr <- paste0("chr",indels$chr)
  res <- tabToIndelsClassification(indels,"testSample","hg38")
  
  expect_count_proportion <- data.frame(sample = "testSample",
                                        del.mh = 7,
                                        del.rep = 4,
                                        del.none = 30,
                                        ins = 8,
                                        complex = 1,
                                        all.del = 41,
                                        all.indels = 50,
                                        del.mh.prop = 0.1707317,
                                        del.rep.prop = 0.09756098,
                                        del.none.prop = 0.7317073,
                                        del.mh.count = 7,
                                        del.rep.count = 4,
                                        del.none.count = 30)
  
  
  
  expect_equal( res$count_proportion, expect_count_proportion,tolerance = 10^-4 )
  
})