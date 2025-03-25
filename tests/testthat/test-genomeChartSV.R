context("testing genome chart SV")

test_that("test genomeChart() using in silico data", {

  filename <- "test.genomeChartSV.pdf"
  resObj <- genomeChartSV(outfilename = filename,
                          SV_bedpe_file = "insilico_sample.bedpe",
                          genome.v = "hg19")
  
  #if no error happen this code can be reached
  expect_true(file.exists(filename))
  
  file.remove(filename)

})

