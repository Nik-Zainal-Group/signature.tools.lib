context("testing Indels classification from VCF file")

test_that("test HRDetect_pipeline() to make sure it raises an error when the input data matrix is not properly formatted.", {

  data_matrix <- NULL
  expect_error(res <- signature.tools.lib::HRDetect_pipeline(data_matrix))
  

  
})

# test_that("test HRDetect_pipeline() runs correctly.", {
#   
#   sample_names <- c("test_hrdetect_1","test_hrdetect_2")
#   col_hrdetect <- c("del.mh.prop", "e.3", "SV3", "SV5", "hrd", "e.8")
#   data_matrix <- matrix(NA,nrow = 2,ncol = 6,dimnames = list(sample_names,col_hrdetect))
#   res <- signature.tools.lib::HRDetect_pipeline(data_matrix)
#   
#   
#   
# })