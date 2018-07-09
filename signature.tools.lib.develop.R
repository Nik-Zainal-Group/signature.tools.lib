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

devtools::document()
setwd("..")
devtools::install("signature.tools.lib")
devtools::test("signature.tools.lib")

