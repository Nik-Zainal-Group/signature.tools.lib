context("testing HRDetect pipeline")

test_that("test HRDetect_pipeline() to make sure it raises an error when the input data matrix is not properly formatted.", {

  data_matrix <- NULL
  expect_null(HRDetect_pipeline(data_matrix))
  

  
})

test_that("test HRDetect_pipeline() runs correctly on one sample. Multiple options set. No bootstrap.", {
  
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
                           optimisation_method = "NNLS",
                           SNV_signature_version = "COSMICv2",
                           bootstrapSignatureFit = FALSE,
                           threshold_percentFit = 10,
                           nparallel = 1)
  #if no error happen this code can be reached
  expect_true(TRUE)
  
})

test_that("test HRDetect_pipeline() runs correctly on two sample. Check random if seed works. No parallel. SigFit bootstrap. No HRDetect bootstrap.", {
  
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
  
  res1 <- HRDetect_pipeline(data_matrix,
                           genome.v = "hg19",
                           SNV_tab_files = SNV_tab_files,
                           SV_bedpe_files = SV_bedpe_files,
                           Indels_vcf_files = Indels_vcf_files,
                           CNV_tab_files = CNV_tab_files,
                           optimisation_method = "NNLS",
                           SNV_signature_version = "COSMICv2",
                           bootstrapSignatureFit = TRUE,
                           nbootFit = 10,
                           threshold_percentFit = 0,
                           nparallel = 1,
                           randomSeed = 1)
  res2 <- HRDetect_pipeline(data_matrix,
                            genome.v = "hg19",
                            SNV_tab_files = SNV_tab_files,
                            SV_bedpe_files = SV_bedpe_files,
                            Indels_vcf_files = Indels_vcf_files,
                            CNV_tab_files = CNV_tab_files,
                            optimisation_method = "NNLS",
                            SNV_signature_version = "COSMICv2",
                            bootstrapSignatureFit = TRUE,
                            nbootFit = 10,
                            threshold_percentFit = 0,
                            nparallel = 1,
                            randomSeed = 1)
  #if no error happen this code can be reached
  expect_equal(res1$hrdetect_output,res2$hrdetect_output)
  
})


test_that("test HRDetect_pipeline() runs correctly on two sample. Check random if seed works. Parallel. SigFit bootstrap. No HRDetect bootstrap.", {
  
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
  
  res1 <- HRDetect_pipeline(data_matrix,
                            genome.v = "hg19",
                            SNV_tab_files = SNV_tab_files,
                            SV_bedpe_files = SV_bedpe_files,
                            Indels_vcf_files = Indels_vcf_files,
                            CNV_tab_files = CNV_tab_files,
                            optimisation_method = "NNLS",
                            organ = "Breast",
                            bootstrapSignatureFit = TRUE,
                            nbootFit = 10,
                            threshold_percentFit = 0,
                            nparallel = 2,
                            randomSeed = 1)
  res2 <- HRDetect_pipeline(data_matrix,
                            genome.v = "hg19",
                            SNV_tab_files = SNV_tab_files,
                            SV_bedpe_files = SV_bedpe_files,
                            Indels_vcf_files = Indels_vcf_files,
                            CNV_tab_files = CNV_tab_files,
                            optimisation_method = "NNLS",
                            organ = "Breast",
                            bootstrapSignatureFit = TRUE,
                            nbootFit = 10,
                            threshold_percentFit = 0,
                            nparallel = 2,
                            randomSeed = 1)
  #if no error happen this code can be reached
  expect_equal(res1$hrdetect_output,res2$hrdetect_output)
  
})

test_that("test HRDetect_pipeline() runs correctly on two sample. Check random if seed works. Parallel. SigFit bootstrap. Use HRDetect bootstrap.", {
  
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
  
  res1 <- HRDetect_pipeline(data_matrix,
                            genome.v = "hg19",
                            SNV_tab_files = SNV_tab_files,
                            SV_bedpe_files = SV_bedpe_files,
                            Indels_vcf_files = Indels_vcf_files,
                            CNV_tab_files = CNV_tab_files,
                            optimisation_method = "NNLS",
                            SNV_signature_version = "COSMICv2",
                            SNV_signature_names = c("Signature1","Signature2","Signature3",
                                                    "Signature5","Signature8","Signature13",
                                                    "Signature17","Signature18"),
                            bootstrapSignatureFit = TRUE,
                            nbootFit = 10,
                            threshold_percentFit = 0,
                            bootstrapHRDetectScores = T,
                            nparallel = 2,
                            randomSeed = 1)
  res2 <- HRDetect_pipeline(data_matrix,
                            genome.v = "hg19",
                            SNV_tab_files = SNV_tab_files,
                            SV_bedpe_files = SV_bedpe_files,
                            Indels_vcf_files = Indels_vcf_files,
                            CNV_tab_files = CNV_tab_files,
                            optimisation_method = "NNLS",
                            SNV_signature_version = "COSMICv2",
                            SNV_signature_names = c("Signature1","Signature2","Signature3",
                                                    "Signature5","Signature8","Signature13",
                                                    "Signature17","Signature18"),
                            bootstrapSignatureFit = TRUE,
                            nbootFit = 10,
                            threshold_percentFit = 0,
                            bootstrapHRDetectScores = T,
                            nparallel = 2,
                            randomSeed = 1)
  #if no error happen this code can be reached
  expect_equal(res1$q_5_50_95,res2$q_5_50_95)
  
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
                           SNV_signature_version = "COSMICv2",
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
                           SNV_signature_version = "COSMICv2",
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
                           organ = "Breast",
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
                           organ = "Breast",
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
                           organ = "Breast",
                           bootstrapHRDetectScores = TRUE,
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
                           organ = "Breast",
                           bootstrapHRDetectScores = TRUE,
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
                           SNV_signature_version = "COSMICv2",
                           SNV_signature_names = c("Signature1","Signature2","Signature3",
                                                   "Signature8","Signature13"),
                           nparallel = 2)
  #if no error happen this code can be reached
  expect_true(TRUE)
  
})

test_that("test HRDetect_pipeline() runs correctly on two samples with subset of COSMIC signatures and TAB indels.", {
  
  sample_names <- c("test_hrdetect_1","test_hrdetect_2")
  col_hrdetect <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
  data_matrix <- matrix(NA,nrow = length(sample_names),ncol = length(col_hrdetect),dimnames = list(sample_names,col_hrdetect))
  
  SNV_tab_files <- c("test_hrdetect_1/test_hrdetect_1.snv.simple.txt",
                     "test_hrdetect_2/test_hrdetect_2.snv.simple.txt")
  names(SNV_tab_files) <- sample_names
  
  SV_bedpe_files <- c("test_hrdetect_1/test_hrdetect_1.sv.bedpe",
                      "test_hrdetect_2/test_hrdetect_2.sv.bedpe")
  names(SV_bedpe_files) <- sample_names
  
  Indels_tab_files <- c("test.tabindels.tsv",
                        "test.tabindels.tsv")
  names(Indels_tab_files) <- sample_names
  
  CNV_tab_files <- c("test_hrdetect_1/test_hrdetect_1.cna.txt",
                     "test_hrdetect_2/test_hrdetect_2.cna.txt")
  names(CNV_tab_files) <- sample_names
  
  res <- HRDetect_pipeline(data_matrix,
                           genome.v = "hg19",
                           SNV_tab_files = SNV_tab_files,
                           SV_bedpe_files = SV_bedpe_files,
                           Indels_tab_files = Indels_tab_files,
                           CNV_tab_files = CNV_tab_files,
                           SNV_signature_version = "COSMICv3.2",
                           SNV_signature_names = c("SBS1","SBS2","SBS3","SBS5",
                                                   "SBS8","SBS13"),
                           nparallel = 2)
  #if no error happen this code can be reached
  expect_true(TRUE)
  
})

test_that("test HRDetect_pipeline() runs correctly on two samples with subset of COSMIC signatures and subset of RefSig Rearr signatures.", {
  
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
                           SNV_signature_version = "COSMICv2",
                           SNV_signature_names = c("Signature1","Signature2","Signature3",
                                                   "Signature8","Signature13"),
                           SV_signature_version = "RefSigv1",
                           SV_signature_names = c("RefSigR1","RefSigR2","RefSigR3",
                                                   "RefSigR4","RefSigR5"),
                           nparallel = 2)
  #if no error happen this code can be reached
  expect_true(TRUE)
  
})

test_that("test HRDetect_pipeline() runs providing signature fit objects from signatureFit_pipeline.", {
  
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
  
  fitPipeline_subs <- signatureFit_pipeline(organ = "Breast",
                                            SNV_tab_files = SNV_tab_files,
                                            useBootstrap = TRUE,
                                            nboot = 100)
  fitPipeline_rearr <- signatureFit_pipeline(organ = "Breast",
                                             SV_bedpe_files = SV_bedpe_files,
                                             useBootstrap = TRUE,
                                             fit_method = "Fit",
                                             nboot = 100)
  
  res <- HRDetect_pipeline(genome.v = "hg19",
                           Indels_vcf_files = Indels_vcf_files,
                           CNV_tab_files = CNV_tab_files,
                           subs_fit_obj = fitPipeline_subs$fitResults,
                           customNameSNV3 = "GEL-Breast_common_SBS3",
                           customNameSNV8 = "GEL-Breast_common_SBS8",
                           rearr_fit_obj = fitPipeline_rearr$fitResults,
                           customNameSV3 = "Breast_G (RefSig R3)",
                           customNameSV5 = "Breast_D (RefSig R5)",
                           organ = "Breast",
                           bootstrapHRDetectScores = TRUE,
                           nparallel = 2)
  #if no error happen this code can be reached
  expect_true(nrow(res$q_5_50_95)==2)
  
})

test_that("test HRDetect_pipeline() runs providing signature fit objects from signatureFit_pipeline and also direct mutations.", {
  
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
  
  fitPipeline_subs <- signatureFit_pipeline(organ = "Breast",
                                            SNV_tab_files = SNV_tab_files[1],
                                            useBootstrap = TRUE,
                                            nboot = 100)
  fitPipeline_rearr <- signatureFit_pipeline(organ = "Breast",
                                             SV_bedpe_files = SV_bedpe_files[1],
                                             useBootstrap = TRUE,
                                             fit_method = "Fit",
                                             nboot = 100)
  
  res <- HRDetect_pipeline(genome.v = "hg19",
                           SNV_tab_files = SNV_tab_files[2],
                           SV_bedpe_files = SV_bedpe_files[2],
                           Indels_vcf_files = Indels_vcf_files,
                           CNV_tab_files = CNV_tab_files,
                           subs_fit_obj = fitPipeline_subs$fitResults,
                           customNameSNV3 = "GEL-Breast_common_SBS3",
                           customNameSNV8 = "GEL-Breast_common_SBS8",
                           rearr_fit_obj = fitPipeline_rearr$fitResults,
                           customNameSV3 = "Breast_G (RefSig R3)",
                           customNameSV5 = "Breast_D (RefSig R5)",
                           organ = "Breast",
                           bootstrapHRDetectScores = TRUE,
                           nparallel = 2)
  #if no error happen this code can be reached
  expect_true(nrow(res$q_5_50_95)==2)
  
})