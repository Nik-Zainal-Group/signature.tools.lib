context("testing conversion from vcf to DNV catalogue")

test_that("test vcfToDNVcatalogue() to convert a VCF file into a DNV catalogue.", {
  
  testDNV_vcffile <- "test.dnv.vcf.bgz"
  
  expected_cat <- read.table(file="test.dnv.vcf.expected.tsv",sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  
  res_cat <- vcfToDNVcatalogue(testDNV_vcffile,genome.v = "hg19")
  obtained_cat <- res_cat$DNV_catalogue
  colnames(obtained_cat) <- "sample"
  
  expect_equal(obtained_cat,expected_cat)
  
})