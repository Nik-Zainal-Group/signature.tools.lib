context("testing Indels classification from data frame")

test_that("test tabToIndelsClassification() using an input file, check count and proportion of indels is correct", {

  indels <- read.table("test.tabindels.tsv",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  
  res <- tabToIndelsClassification(indels,"testSample","hg19")
  
  expect_count_proportion <- data.frame(sample = "testSample",
                                        del.mh = 22,
                                        del.rep = 15,
                                        del.none = 4,
                                        ins = 8,
                                        all.indels = 41,
                                        del.mh.prop = 0.536585365853659,
                                        del.rep.prop = 0.365853658536585,
                                        del.none.prop = 0.0975609756097561,
                                        del.mh.count = 22,
                                        del.rep.count = 15,
                                        del.none.count = 4)

  
  
  expect_equal( res$count_proportion, expect_count_proportion )
  
})