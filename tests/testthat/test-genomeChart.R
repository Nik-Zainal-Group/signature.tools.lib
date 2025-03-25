context("testing genome chart with circlize on subs, indels, CNVs and rearr data")

test_that("test genomeChart() using subs tab, indels vcfs, copy number, and rearrangment bedpe input file formats", {

  resObj <- genomeChart(outfilename = "test.genomeChart.pdf",
                        SNV_tab_file = "test_hrdetect_1/test_hrdetect_1.snv.simple.txt",
                        Indels_vcf_file = "test_hrdetect_1/test_hrdetect_1.indel.vcf.gz",
                        CNV_tab_file = "test_hrdetect_1/test_hrdetect_1.cna.txt",
                        SV_bedpe_file = "test_hrdetect_1/test_hrdetect_1.sv.bedpe",
                        sample_name = "test")

  file.remove("test.genomeChart.pdf")

  #if no error happen this code can be reached
  expect_true(TRUE)

})

