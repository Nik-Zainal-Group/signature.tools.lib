context("testing conversion from dataframe of SNV to SNV catalogue")

test_that("test vcfToSNVcatalogue() to convert a VCF file with SNV into an SNV catalogue.", {
  
  testSNV_vcffile <- "test.sub.vcf.gz"
  
  expected_cat <- read.table(file="test2.snv.tab",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  
  res_cat <- vcfToSNVcatalogue(testSNV_vcffile,genome.v = "hg19")
  
  expect_equal(res_cat$catalogue,expected_cat)
  
  #expect_error(res <- HRDetect_pipeline(data_matrix))
  #write.table(res_cat$catalogue,file="test2.snv.tab",sep = "\t",col.names = TRUE,row.names = TRUE,quote = FALSE)
  
})