context("testing Indels classification from VCF file")

test_that("test vcfToIndelsClassification() using an input file, check count and proportion of indels is correct", {

  res <- vcfToIndelsClassification("test.indel.vcf.gz","testSample","hg19")
  
  expect_count_proportion <- data.frame(sample = "testSample",
                                        del.mh = 48,
                                        del.rep = 193,
                                        del.none = 70,
                                        ins = 146,
                                        complex = 3,
                                        all.del = 311,
                                        all.indels = 460,
                                        del.mh.prop = 0.1543408360,
                                        del.rep.prop = 0.6205787781,
                                        del.none.prop = 0.2250803859,
                                        del.mh.count = 48,
                                        del.rep.count = 193,
                                        del.none.count = 70)

  
  
  expect_equal( res$count_proportion, expect_count_proportion )
  
})

test_that("test vcfToIndelsClassification() using an input file, check count and proportion of indels is correct", {
  
  res <- vcfToIndelsClassification("test.indel.vcf.gz","testSample","hg38")
  
  expect_count_proportion <- data.frame(sample = "testSample",
                                        del.mh = 20,
                                        del.rep = 67,
                                        del.none = 224,
                                        ins = 146,
                                        complex = 3,
                                        all.del = 311,
                                        all.indels = 460,
                                        del.mh.prop = 20/311,
                                        del.rep.prop = 67/311,
                                        del.none.prop = 224/311,
                                        del.mh.count = 20,
                                        del.rep.count = 67,
                                        del.none.count = 224)
  
  
  
  expect_equal( res$count_proportion, expect_count_proportion )
  
})