context("testing building a rearrangement catalogue from a BEDPE file")

test_that("test bedpeToRearrCatalogue() using a simple input file 1/2", {

  sv_bedpe <- read.table("test.bedpe",sep = "\t",header = TRUE,
                         stringsAsFactors = FALSE,check.names = FALSE)
  sv_bedpe <- sv_bedpe[,c("chrom1","start1","end1","chrom2","start2","end2","sample","strand1","strand2")]
  expected_res <- read.table("test.cat",sep = "\t",header = TRUE,
                             stringsAsFactors = FALSE,check.names = FALSE)
  reslist <- bedpeToRearrCatalogue(sv_bedpe)
  res <- reslist$rearr_catalogue
  
  expect_equal( res, expected_res )

})

test_that("test bedpeToRearrCatalogue() using a simple input file 2/2", {
  
  sv_bedpe <- read.table("test_hrdetect_1/test_hrdetect_1.sv.bedpe",sep = "\t",header = TRUE,
                         stringsAsFactors = FALSE,check.names = FALSE)
  sv_bedpe <- sv_bedpe[,c("chrom1","start1","end1","chrom2","start2","end2","strand1","strand2")]
  sv_bedpe$sample <- "Freq"
  expected_res <- read.table("test_hrdetect_1/test_hrdetect_1_expected_rearrcat.txt",sep = "\t",header = TRUE,
                             stringsAsFactors = FALSE,check.names = FALSE)
  reslist <- bedpeToRearrCatalogue(sv_bedpe)
  res <- reslist$rearr_catalogue
  
  expect_equal( res, expected_res )
  
})