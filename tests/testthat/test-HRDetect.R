context("testing Indels classification from VCF file")

test_that("test HRDetect_pipeline() to make sure it raises an error when the input data matrix is not properly formatted.", {

  data_matrix <- NULL
  expect_error(res <- signature.tools.lib::HRDetect_pipeline(data_matrix))
  

  
})

test_that("test HRDetect_pipeline() runs correctly.", {

  sample_names <- c("test_hrdetect_1","test_hrdetect_2")
  col_hrdetect <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
  data_matrix <- matrix(NA,nrow = 2,ncol = 6,dimnames = list(sample_names,col_hrdetect))
  
  SNV_tab_files <- c("test_hrdetect_1/test_hrdetect_1.snv.simple.txt",
    "test_hrdetect_2/test_hrdetect_2.snv.simple.txt")
  names(SNV_tab_files) <- sample_names
  
  res <- signature.tools.lib::HRDetect_pipeline(data_matrix,
                                                genome.v="hg19",
                                                SNV_tab_files=SNV_tab_files,
                                                nparallel = 1)
  res$data_matrix


})