#file with commands to set up the signature-tools-lib R package
#setwd("~/sandbox/git")

# install.packages("devtools")
# install.packages("roxygen2")
# 
# devtools::create("signature.tools.lib")

#setwd("~/sandbox/git/signature.tools.lib")

devtools::use_package("VariantAnnotation")
devtools::use_package("SummarizedExperiment")
devtools::use_package("BiocGenerics")
devtools::use_package("GenomeInfoDb")
devtools::use_package("BSgenome")
devtools::use_package("BSgenome.Hsapiens.UCSC.hg38")
devtools::use_package("BSgenome.Hsapiens.1000genomes.hs37d5")
devtools::use_package("BSgenome.Mmusculus.UCSC.mm10")
devtools::use_package("BSgenome.Cfamiliaris.UCSC.canFam3")
devtools::use_package("NMF")
devtools::use_package("foreach")
devtools::use_package("doParallel")
devtools::use_package("doMC")
devtools::use_package("lpSolve")
devtools::use_package("ggplot2")
devtools::use_package("methods")
devtools::use_package("cluster")
devtools::use_package("stats")
devtools::use_package("NNLM")
devtools::use_package("nnls")
devtools::use_package("GenSA")
devtools::use_package("gmp")
devtools::use_package("plyr")
devtools::use_package("RCircos")
devtools::use_package("scales")
devtools::use_package("GenomicRanges")
devtools::use_package("IRanges")
devtools::load_all()

#add internal data
RS.Breast560 <- read.table("data/Breast560_rearrangement.signatures.txt", sep="\t", header=T, as.is=T, check.names = FALSE)
cosmic30 <- read.table("data/COSMIC30_subs_signatures.txt", sep="\t", header=T, as.is=T, check.names = FALSE)
all_organ_sigs_subs <- read.table("data/2019_01_10_all_PCAWG_sigs_subs.tsv", sep="\t", header=T, as.is=T, check.names = FALSE)
all_organ_sigs_rearr <- read.table("data/2019_01_10_all_PCAWG_sigs_rearr.tsv", sep="\t", header=T, as.is=T, check.names = FALSE)
conversion_matrix_subs <- read.table("data/2019_01_10_ConversionMatrix_subs.tsv", sep="\t", header=T, as.is=T, check.names = FALSE)
conversion_matrix_rearr <- read.table("data/2019_01_10_ConversionMatrix_rearr.tsv", sep="\t", header=T, as.is=T, check.names = FALSE)
load("data/chrominfo.RData")
load("data/chrominfo.snp6.hg19.RData")
usethis::use_data(RS.Breast560,
                   cosmic30,
                   chrominfo,
                   chrominfo.snp6,
                   all_organ_sigs_subs,
                   all_organ_sigs_rearr,
                   conversion_matrix_subs,
                   conversion_matrix_rearr,
                   internal = TRUE,overwrite = TRUE)

devtools::document()
devtools::install()

#test all
devtools::test()

#some individual tests
devtools::test(pkg = ".",filter = "ascatToHRDLOH")
devtools::test(pkg = ".",filter = "tabToIndelsClassification")
devtools::test(pkg = ".",filter = "vcfToIndelsClassification")
devtools::test(pkg = ".",filter = "HRDetect")
devtools::test(pkg = ".",filter = "tabToSNVcatalogue")
devtools::test(pkg = ".",filter = "vcfToSNVcatalogue")
devtools::test(pkg = ".",filter = "bedpeToRearrCatalogue")
devtools::test(pkg = ".",filter = "genomePlot")
devtools::test(pkg = ".",filter = "SignatureExtraction")
devtools::test(pkg = ".",filter = "SignatureFit")



