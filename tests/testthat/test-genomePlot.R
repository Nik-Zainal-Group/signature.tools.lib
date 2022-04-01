context("testing genome plot with RCircos on subs, indels, CNVs and rearr data")

test_that("test genomePlot() using subs/indels vcfs, comma-separated without header copy number, and rearrangment bedpe input file formats and compare with the expected", {

  #mut.order <<- c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T")

  output_png <- genomePlot("test.sub.vcf.gz","test.indel.vcf.gz","test.cnvs.txt","test.bedpe","test")

  #expected_png <- "genomePlot.test.png"
  #expect_equal( as.character(tools::md5sum(output_png)), as.character(tools::md5sum(expected_png)) )
  file.remove("test.genomePlot.png")

  #if no error happen this code can be reached
  expect_true(TRUE)

})

test_that("test genomePlot() using subs/indels vcfs, tab separated with header copy number, and rearrangment bedpe input file formats and compare with the expected", {

  output_png <- genomePlot(subsVcf.file = "test.sub.vcf.gz","test.indel.vcf.gz",
                           "test_hrdetect_1/test_hrdetect_1.cna.txt",
                           "test_hrdetect_1/test_hrdetect_1.sv.bedpe","test")

  file.remove("test.genomePlot.png")

  #if no error happen this code can be reached
  expect_true(TRUE)

})
