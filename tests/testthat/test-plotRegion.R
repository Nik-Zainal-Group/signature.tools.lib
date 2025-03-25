context("testing plotRegion on subs, indels, CNVs and rearr data")

test_that("test plotRegion() using subs tab, indels vcfs, copy number, and rearrangment bedpe input file formats", {

  outfile <- "plotRegion_test.pdf"
  res_plot <- plotRegion(outfilename = outfile,
                         sample_name = "test sample",
                         chr = "1",
                         rstart = 1,
                         rend = 240000000,
                         plotSNVandIndels = TRUE,
                         SNV_vcf_file = NULL,
                         SNV_tab_file = "test_hrdetect_2/test_hrdetect_2.snv.simple.txt",
                         snvs_table = NULL,
                         Indels_vcf_file = "test_hrdetect_2/test_hrdetect_2.indel.vcf.gz",
                         Indels_tab_file = NULL,
                         indels_table = NULL,
                         CNV_tab_file = "test_hrdetect_2/test_hrdetect_2.cna.txt",
                         CNV_table = NULL,
                         SV_bedpe_file = "test_hrdetect_2/test_hrdetect_2.sv.bedpe",
                         sv_table = NULL,
                         plot_title = "test sample",
                         spreadSVlabels = TRUE,
                         genome.v = "hg19",
                         debug=FALSE)

  file.remove(outfile)

  #if no error happen this code can be reached
  expect_true(TRUE)

})


test_that("test plotRegion() using copy number, and rearrangment bedpe input file formats", {
  
  outfile <- "plotRegion_test_SNV_SV.pdf"
  res_plot <- plotRegion(outfilename = outfile,
                         sample_name = "test sample",
                         chr = "1",
                         rstart = 1,
                         rend = 240000000,
                         plotSNVandIndels = FALSE,
                         CNV_tab_file = "test_hrdetect_2/test_hrdetect_2.cna.txt",
                         CNV_table = NULL,
                         SV_bedpe_file = "test_hrdetect_2/test_hrdetect_2.sv.bedpe",
                         sv_table = NULL,
                         plot_title = "test sample",
                         spreadSVlabels = TRUE,
                         genome.v = "hg19",
                         debug=FALSE)
  
  file.remove(outfile)
  
  #if no error happen this code can be reached
  expect_true(TRUE)
  
})

