context("testing HRDetect pipeline")

test_that("test HRDetect_pipeline() to make sure it raises an error when the input data matrix is not properly formatted.", {

  data_matrix <- NULL
  expect_error(res <- HRDetect_pipeline(data_matrix))
  

  
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
  
  Indels_vcf_files <- c("test_hrdetect_1/test_hrdetect_1.indel.vcf.gz",
                      "test_hrdetect_2/test_hrdetect_2.indel.vcf.gz")
  names(Indels_vcf_files) <- sample_names
  
  CNV_tab_files <- c("test_hrdetect_1/test_hrdetect_1.cna.txt",
                     "test_hrdetect_2/test_hrdetect_2.cna.txt")
  names(CNV_tab_files) <- sample_names
  
  res <- HRDetect_pipeline(data_matrix,
                                                genome.v = "hg19",
                                                SNV_tab_files = SNV_tab_files,
                                                SV_bedpe_files = SV_bedpe_files,
                                                Indels_vcf_files = Indels_vcf_files,
                                                CNV_tab_files = CNV_tab_files,
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
  
  Indels_vcf_files <- c("test_hrdetect_1/test_hrdetect_1.indel.vcf.gz")
  names(Indels_vcf_files) <- sample_names
  
  CNV_tab_files <- c("test_hrdetect_1/test_hrdetect_1.cna.txt")
  names(CNV_tab_files) <- sample_names
  
  res <- HRDetect_pipeline(data_matrix,
                                                genome.v = "hg19",
                                                SNV_tab_files = SNV_tab_files,
                                                SV_bedpe_files = SV_bedpe_files,
                                                Indels_vcf_files = Indels_vcf_files,
                                                CNV_tab_files = CNV_tab_files,
                                                nparallel = 1)
  #if no error happen this code can be reached
  expect_true(TRUE)
  
})

test_that("test HRDetect_pipeline() runs correctly on one sample with Breast signatures and no data_matrix provided.", {
  
  sample_names <- c("test_hrdetect_1")
  
  SNV_tab_files <- c("test_hrdetect_1/test_hrdetect_1.snv.simple.txt")
  names(SNV_tab_files) <- sample_names
  
  SV_bedpe_files <- c("test_hrdetect_1/test_hrdetect_1.sv.bedpe")
  names(SV_bedpe_files) <- sample_names
  
  Indels_vcf_files <- c("test_hrdetect_1/test_hrdetect_1.indel.vcf.gz")
  names(Indels_vcf_files) <- sample_names
  
  CNV_tab_files <- c("test_hrdetect_1/test_hrdetect_1.cna.txt")
  names(CNV_tab_files) <- sample_names
  
  res <- HRDetect_pipeline(genome.v = "hg19",
                           SNV_tab_files = SNV_tab_files,
                           SV_bedpe_files = SV_bedpe_files,
                           Indels_vcf_files = Indels_vcf_files,
                           CNV_tab_files = CNV_tab_files,
                           signature_type = "Breast",
                           nparallel = 1)
  #if no error happen this code can be reached
  expect_true(TRUE)
  
})


test_that("test HRDetect_pipeline() runs correctly on two samples with Breast signatures and no data_matrix provided.", {
  
  sample_names <- c("test_hrdetect_1","test_hrdetect_2")
  
  SNV_tab_files <- c("test_hrdetect_1/test_hrdetect_1.snv.simple.txt",
                     "test_hrdetect_2/test_hrdetect_2.snv.simple.txt")
  names(SNV_tab_files) <- sample_names
  
  SV_bedpe_files <- c("test_hrdetect_1/test_hrdetect_1.sv.bedpe",
                      "test_hrdetect_2/test_hrdetect_2.sv.bedpe")
  names(SV_bedpe_files) <- sample_names
  
  Indels_vcf_files <- c("test_hrdetect_1/test_hrdetect_1.indel.vcf.gz",
                        "test_hrdetect_2/test_hrdetect_2.indel.vcf.gz")
  names(Indels_vcf_files) <- sample_names
  
  CNV_tab_files <- c("test_hrdetect_1/test_hrdetect_1.cna.txt",
                     "test_hrdetect_2/test_hrdetect_2.cna.txt")
  names(CNV_tab_files) <- sample_names
  
  res <- HRDetect_pipeline(genome.v = "hg19",
                           SNV_tab_files = SNV_tab_files,
                           SV_bedpe_files = SV_bedpe_files,
                           Indels_vcf_files = Indels_vcf_files,
                           CNV_tab_files = CNV_tab_files,
                           signature_type = "Breast",
                           nparallel = 2)
  #if no error happen this code can be reached
  expect_true(TRUE)
  
})

test_that("test HRDetect_pipeline() runs correctly on one sample with Breast signatures and no data_matrix provided and with bootstrap.", {
  
  sample_names <- c("test_hrdetect_1")
  
  SNV_tab_files <- c("test_hrdetect_1/test_hrdetect_1.snv.simple.txt")
  names(SNV_tab_files) <- sample_names
  
  SV_bedpe_files <- c("test_hrdetect_1/test_hrdetect_1.sv.bedpe")
  names(SV_bedpe_files) <- sample_names
  
  Indels_vcf_files <- c("test_hrdetect_1/test_hrdetect_1.indel.vcf.gz")
  names(Indels_vcf_files) <- sample_names
  
  CNV_tab_files <- c("test_hrdetect_1/test_hrdetect_1.cna.txt")
  names(CNV_tab_files) <- sample_names
  
  res <- HRDetect_pipeline(genome.v = "hg19",
                           SNV_tab_files = SNV_tab_files,
                           SV_bedpe_files = SV_bedpe_files,
                           Indels_vcf_files = Indels_vcf_files,
                           CNV_tab_files = CNV_tab_files,
                           signature_type = "Breast",
                           bootstrap_scores = TRUE,
                           nparallel = 1)
  
  #if no error happen this code can be reached
  expect_true(nrow(res$q_5_50_95)==1)
  
})

test_that("test HRDetect_pipeline() runs correctly on two samples with Breast signatures and no data_matrix provided and with bootstrap.", {
  
  sample_names <- c("test_hrdetect_1","test_hrdetect_2")
  
  SNV_tab_files <- c("test_hrdetect_1/test_hrdetect_1.snv.simple.txt",
                     "test_hrdetect_2/test_hrdetect_2.snv.simple.txt")
  names(SNV_tab_files) <- sample_names
  
  SV_bedpe_files <- c("test_hrdetect_1/test_hrdetect_1.sv.bedpe",
                      "test_hrdetect_2/test_hrdetect_2.sv.bedpe")
  names(SV_bedpe_files) <- sample_names
  
  Indels_vcf_files <- c("test_hrdetect_1/test_hrdetect_1.indel.vcf.gz",
                        "test_hrdetect_2/test_hrdetect_2.indel.vcf.gz")
  names(Indels_vcf_files) <- sample_names
  
  CNV_tab_files <- c("test_hrdetect_1/test_hrdetect_1.cna.txt",
                     "test_hrdetect_2/test_hrdetect_2.cna.txt")
  names(CNV_tab_files) <- sample_names
  
  res <- HRDetect_pipeline(genome.v = "hg19",
                           SNV_tab_files = SNV_tab_files,
                           SV_bedpe_files = SV_bedpe_files,
                           Indels_vcf_files = Indels_vcf_files,
                           CNV_tab_files = CNV_tab_files,
                           signature_type = "Breast",
                           bootstrap_scores = TRUE,
                           nparallel = 2)
  #if no error happen this code can be reached
  expect_true(nrow(res$q_5_50_95)==2)
  
})

test_that("test HRDetect_pipeline() runs correctly on two samples with subset of COSMIC signatures.", {
  
  sample_names <- c("test_hrdetect_1","test_hrdetect_2")
  col_hrdetect <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
  data_matrix <- matrix(NA,nrow = length(sample_names),ncol = length(col_hrdetect),dimnames = list(sample_names,col_hrdetect))
  
  SNV_tab_files <- c("test_hrdetect_1/test_hrdetect_1.snv.simple.txt",
                     "test_hrdetect_2/test_hrdetect_2.snv.simple.txt")
  names(SNV_tab_files) <- sample_names
  
  SV_bedpe_files <- c("test_hrdetect_1/test_hrdetect_1.sv.bedpe",
                      "test_hrdetect_2/test_hrdetect_2.sv.bedpe")
  names(SV_bedpe_files) <- sample_names
  
  Indels_vcf_files <- c("test_hrdetect_1/test_hrdetect_1.indel.vcf.gz",
                        "test_hrdetect_2/test_hrdetect_2.indel.vcf.gz")
  names(Indels_vcf_files) <- sample_names
  
  CNV_tab_files <- c("test_hrdetect_1/test_hrdetect_1.cna.txt",
                     "test_hrdetect_2/test_hrdetect_2.cna.txt")
  names(CNV_tab_files) <- sample_names
  
  res <- HRDetect_pipeline(data_matrix,
                           genome.v = "hg19",
                           SNV_tab_files = SNV_tab_files,
                           SV_bedpe_files = SV_bedpe_files,
                           Indels_vcf_files = Indels_vcf_files,
                           CNV_tab_files = CNV_tab_files,
                           signature_type = "COSMIC",
                           cosmic_siglist = c(1,2,3,8,13),
                           nparallel = 2)
  #if no error happen this code can be reached
  expect_true(TRUE)
  
})
