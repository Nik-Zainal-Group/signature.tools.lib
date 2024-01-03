context("testing assign signature probability to mutations")

test_that("test assignSignatureProbabilityToMutations() using SNVs.", {
  
  expectedResult <- readTable("../testthat/assignSigProbability/sigsProb_test_hrdetect_1.tsv")
  
  testSNVdf <- read.table(file = "test_hrdetect_1/test_hrdetect_1.snv.simple.txt",
                          sep = "\t",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
  
  snvres <- tabToSNVcatalogue(testSNVdf,genome.v = "hg19")
  
  signames <- c("SBS1","SBS2","SBS3","SBS5","SBS8","SBS18")
  signatures <- getReferenceSignatures()
  signatures <- signatures[,signames,drop=F]
  
  sampleMutations <- snvres$muts
  
  # resfit <- Fit(catalogues = snvres$catalogue,signatures = signatures)
  # 
  # sampleSigsExposures <- resfit$exposures
  sampleSigsExposures <- matrix(c(819.2134,178.9108,255.4519,1107.957,216.8609,392.606),nrow = 1, ncol = ncol(signatures),dimnames = list("sample1",colnames(signatures)))
  
  reassign <- assignSignatureProbabilityToMutations(sampleMutations = sampleMutations,
                                                    sampleSigsExposures = sampleSigsExposures,
                                                    signatures = signatures)
  
  expect_equal(reassign,expectedResult)
  
})