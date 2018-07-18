context("testing computation of HRD-LOH from ascat file")

test_that("test ascatToHRDLOH() using an input file 1/2", {
  
  ascat_test.df <- read.table("ascat_test.csv", sep=",", header=F, as.is=T, check.names = FALSE)
  colnames(ascat_test.df) <- c('seg_no', 'Chromosome', 'chromStart', 'chromEnd', 'total.copy.number.inNormal', 'minor.copy.number.inNormal', 'total.copy.number.inTumour', 'minor.copy.number.inTumour')
  
  res <- ascatToHRDLOH(ascat_test.df,"testSample")
  
  expected_res <- 19
  names(expected_res) <- "testSample"
  
  expect_equal( res, expected_res )
  
})

test_that("test ascatToHRDLOH() using an input file 2/2", {
  
  ascat_test.df <- read.table("test_hrdetect_1/test_hrdetect_1.cna.txt", sep="\t", header=T, as.is=T, check.names = FALSE)

  #old code to convert from pan can to necessary input format
  # ascat_test.df <- ascat_test.df[complete.cases(ascat_test.df),]
  # ascat_test.df$seg_no <- 1:nrow(ascat_test.df)
  # ascat_test.df$Chromosome <- ascat_test.df$chromosome
  # ascat_test.df$chromStart <- ascat_test.df$start
  # ascat_test.df$chromEnd <- ascat_test.df$end
  # ascat_test.df$total.copy.number.inNormal <- 2
  # ascat_test.df$minor.copy.number.inNormal <- 1
  # ascat_test.df$total.copy.number.inTumour <- ascat_test.df$total_cn
  # ascat_test.df$minor.copy.number.inTumour <- ascat_test.df$minor_cn
  # 
  # ascat_test.df <- ascat_test.df[,c('seg_no', 'Chromosome', 'chromStart', 'chromEnd', 'total.copy.number.inNormal', 'minor.copy.number.inNormal', 'total.copy.number.inTumour', 'minor.copy.number.inTumour')]
  # 
  # write.table(ascat_test.df,file = "test_hrdetect_1/test_hrdetect_1.cna.txt", sep="\t",quote = FALSE,row.names = FALSE,col.names = TRUE)
  # 
  res <- ascatToHRDLOH(ascat_test.df,"test_hrdetect_1")
  
  expected_res <- 5
  names(expected_res) <- "test_hrdetect_1"
  
  expect_equal( res, expected_res )
  
})