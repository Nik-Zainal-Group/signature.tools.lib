context("testing HRDetect pipeline")

test_that("test HRDetect_pipeline() to make sure it raises an error when the input data matrix is not properly formatted.", {

  data_matrix <- NULL
  expect_error(res <- signature.tools.lib::HRDetect_pipeline(data_matrix))
  

  
})

test_that("test HRDetect_pipeline() runs correctly on two samples.", {

  sample_names <- c("test_hrdetect_1","test_hrdetect_2")
  col_hrdetect <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
  data_matrix <- matrix(NA,nrow = length(sample_names),ncol = length(col_hrdetect),dimnames = list(sample_names,col_hrdetect))
  
  SNV_tab_files <- c("test_hrdetect_1/test_hrdetect_1.snv.simple.txt",
                     "test_hrdetect_2/test_hrdetect_2.snv.simple.txt")
  names(SNV_tab_files) <- sample_names
  
  SV_bedpe_files <- c("test_hrdetect_1/test_hrdetect_1.sv.bedpe",
                      "test_hrdetect_2/test_hrdetect_2.sv.bedpe")
  names(SV_bedpe_files) <- sample_names
  
  res <- signature.tools.lib::HRDetect_pipeline(data_matrix,
                                                genome.v="hg19",
                                                SNV_tab_files=SNV_tab_files,
                                                SV_bedpe_files=SV_bedpe_files,
                                                nparallel = 2)
  #if no error happen this code can be reached
  expect_true(TRUE)

})

test_that("test HRDetect_pipeline() runs correctly on one sample.", {
  
  sample_names <- c("test_hrdetect_1")
  col_hrdetect <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
  data_matrix <- matrix(NA,nrow = length(sample_names),ncol = length(col_hrdetect),dimnames = list(sample_names,col_hrdetect))
  
  SNV_tab_files <- c("test_hrdetect_1/test_hrdetect_1.snv.simple.txt")
  names(SNV_tab_files) <- sample_names
  
  SV_bedpe_files <- c("test_hrdetect_1/test_hrdetect_1.sv.bedpe")
  names(SV_bedpe_files) <- sample_names
  
  res <- signature.tools.lib::HRDetect_pipeline(data_matrix,
                                                genome.v="hg19",
                                                SNV_tab_files=SNV_tab_files,
                                                SV_bedpe_files=SV_bedpe_files,
                                                nparallel = 1)
  #if no error happen this code can be reached
  expect_true(TRUE)
  
})
