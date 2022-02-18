context("testing SignatureFit pipeline")

test_that("test that signatureFit_pipeline() works on single substitutions with default parameters", {
  
  snv_files <- c("test_hrdetect_1/test_hrdetect_1.snv.simple.txt",
                 "test_hrdetect_2/test_hrdetect_2.snv.simple.txt")
  names(snv_files) <- c("sample1","sample2")
  res <- signatureFit_pipeline(SNV_tab_files = snv_files,
                                       organ = "Breast")
  
  #here goes the test
  expect_true(!is.null(res$fitResults))
})

test_that("test that signatureFit_pipeline() works on double substitutions with default parameters", {
  
  snv_files <- c("test_hrdetect_1/test_hrdetect_1.snv.simple.txt",
                 "test_hrdetect_2/test_hrdetect_2.snv.simple.txt")
  names(snv_files) <- c("sample1","sample2")
  res <- signatureFit_pipeline(DNV_tab_files = snv_files,
                               fit_method = "Fit",
                               organ = "Breast")
  
  #here goes the test
  expect_true(!is.null(res$fitResults))
})

test_that("test that signatureFit_pipeline() works on rearrangements with default parameters", {
  
  sv_files <- c("test_hrdetect_1/test_hrdetect_1.sv.bedpe",
                 "test_hrdetect_2/test_hrdetect_2.sv.bedpe")
  names(sv_files) <- c("sample1","sample2")
  res <- signatureFit_pipeline(SV_bedpe_files = sv_files,
                               fit_method = "Fit",
                               organ = "Breast")
  
  #here goes the test
  expect_true(!is.null(res$fitResults))
})

test_that("test that signatureFit_pipeline() works on single substitutions with all COSMICv2 signatures", {
  
  snv_files <- c("test_hrdetect_1/test_hrdetect_1.snv.simple.txt",
                 "test_hrdetect_2/test_hrdetect_2.snv.simple.txt")
  names(snv_files) <- c("sample1","sample2")
  res <- signatureFit_pipeline(SNV_tab_files = snv_files,
                               fit_method = "Fit",
                               signature_version = "COSMICv2")
  
  #here goes the test
  expect_true(!is.null(res$fitResults))
})

test_that("test that signatureFit_pipeline() works on single substitutions with selected COSMICv2 signatures", {
  
  snv_files <- c("test_hrdetect_1/test_hrdetect_1.snv.simple.txt",
                 "test_hrdetect_2/test_hrdetect_2.snv.simple.txt")
  names(snv_files) <- c("sample1","sample2")
  res <- signatureFit_pipeline(SNV_tab_files = snv_files,
                               fit_method = "Fit",
                               signature_version = "COSMICv2",
                               signature_names = paste0("Signature",c(1,2,3,5,8,13,17,18)))
  
  #here goes the test
  expect_true(!is.null(res$fitResults))
})

test_that("test that signatureFit_pipeline() works on single substitutions with selected COSMICv3.2 signatures", {
  
  snv_files <- c("test_hrdetect_1/test_hrdetect_1.snv.simple.txt",
                 "test_hrdetect_2/test_hrdetect_2.snv.simple.txt")
  names(snv_files) <- c("sample1","sample2")
  res <- signatureFit_pipeline(SNV_tab_files = snv_files,
                               fit_method = "Fit",
                               signature_version = "COSMICv3.2",
                               signature_names = paste0("SBS",c("1","2","3","5","8","13","17a","17b","18")))
  
  #here goes the test
  expect_true(!is.null(res$fitResults))
})
