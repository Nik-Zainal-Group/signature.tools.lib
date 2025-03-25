context("testing genome chart with circlize")

test_that("test genomeChart() using subs tab, indels vcfs, copy number, and rearrangment bedpe input file formats", {

  filename <- "test.genomeChart.pdf"
  resObj <- genomeChart(outfilename = filename,
                        SNV_tab_file = "test_hrdetect_1/test_hrdetect_1.snv.simple.txt",
                        Indels_vcf_file = "test_hrdetect_1/test_hrdetect_1.indel.vcf.gz",
                        CNV_tab_file = "test_hrdetect_1/test_hrdetect_1.cna.txt",
                        SV_bedpe_file = "test_hrdetect_1/test_hrdetect_1.sv.bedpe",
                        sample_name = "test")
  
  #if no error happen this code can be reached
  expect_true(file.exists(filename))
  
  file.remove(filename)

})

