context("testing Indels classification from VCF file")

test_that("test vcfToIndelsClassification() using an input file, check count and proportion of indels is correct", {

  res <- vcfToIndelsClassification("test.indel.vcf.gz","testSample","hg19")
  
  expect_count_proportion <- data.frame(sample = "testSample",
                                        del.mh = 48,
                                        del.rep = 193,
                                        del.none = 70,
                                        ins = 146,
                                        all.indels = 311,
                                        del.mh.prop = 0.1543408360,
                                        del.rep.prop = 0.6205787781,
                                        del.none.prop = 0.2250803859,
                                        del.mh.count = 48,
                                        del.rep.count = 193,
                                        del.none.count = 70)

  
  
  expect_equal( res$count_proportion, expect_count_proportion )
  
})