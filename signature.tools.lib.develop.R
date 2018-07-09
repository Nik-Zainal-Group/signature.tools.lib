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
devtools::use_package("BSgenome.Hsapiens.NCBI.GRCh38")
devtools::use_package("BSgenome.Hsapiens.UCSC.hg19")
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

#add data
RS.Breast560 <- read.table("data/Breast560_rearrangement.signatures.txt", sep="\t", header=T, as.is=T, check.names = FALSE)
cosmic30 <- read.table("data/COSMIC30_subs_signatures.txt", sep="\t", header=T, as.is=T, check.names = FALSE)
devtools::use_data(RS.Breast560,cosmic30, internal = TRUE,overwrite = TRUE)

devtools::document()
devtools::install()
devtools::test()

