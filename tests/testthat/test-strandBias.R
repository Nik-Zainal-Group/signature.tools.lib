context("testing strand bias functions")

test_that("test sampleStrandBias()", {
  
  expectedResult <- readTable("../testthat/strandBias/mmr1_biasTable.tsv")
  
  mmr1 <- readTable("../testthat/strandBias/mmr_muts_1.tsv")
  mmr1_sb <- sampleStrandBias(snv_table = mmr1,genomev = "hg19")
  
  expect_equal(mmr1_sb$bias_results_single,expectedResult)
  
})

test_that("test sampleStrandBias(), combineStrandBiasResults() and plotStrandBiasResults()", {
  
  mmr1 <- readTable("../testthat/strandBias/mmr_muts_1.tsv")
  mmr1_sb <- sampleStrandBias(snv_table = mmr1,genomev = "hg19")
  mmr2 <- readTable("../testthat/strandBias/mmr_muts_2.tsv")
  mmr2_sb <- sampleStrandBias(snv_table = mmr2,genomev = "hg19")
  mmr3 <- readTable("../testthat/strandBias/mmr_muts_3.tsv")
  mmr3_sb <- sampleStrandBias(snv_table = mmr3,genomev = "hg19")

  # need to try to combine also
  biasResObjList <- list()
  biasResObjList[["MMRd sample 1"]] <- mmr1_sb
  biasResObjList[["MMRd sample 2"]] <- mmr2_sb
  biasResObjList[["MMRd sample 3"]] <- mmr3_sb
  
  res_CombMMR <- combineStrandBiasResults(biasResObjList = biasResObjList)
  plotStrandBiasResults(biasResObj = res_CombMMR,
                        filename = "combinedSamples_test_plot_MMR.pdf",
                        addToTitle = ", SBS44 samples")
  
  # clean up
  unlink("combinedSamples_test_plot_MMR.pdf")
  unlink("combinedSamples_test_plot_MMR_ratio_replication_bias_meanStErr.pdf")
  unlink("combinedSamples_test_plot_MMR_ratio_transcription_bias_meanStErr.pdf")
  
  #here goes the test
  expect_true(TRUE)
  
})