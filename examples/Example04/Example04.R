#set directory to this file location

#import the package
library(signature.tools.lib)

#set sample names
sample_names <- c("sample1","sample2")

#set the file names.
SNV_tab_files <- c("../../tests/testthat/test_hrdetect_1/test_hrdetect_1.snv.simple.txt",
                   "../../tests/testthat/test_hrdetect_2/test_hrdetect_2.snv.simple.txt")

#name the vectors entries with the sample names
names(SNV_tab_files) <- sample_names

#use the signatureFit_pipeline function to generate catalogues and perform
#signature fit all at once
pipeline_subs_res <- signatureFit_pipeline(SNV_tab_files = SNV_tab_files,
                                           organ = "Breast",genome.v = "hg19",
                                           fit_method = "FitMS",nparallel = 2)
#plot the catalogues. The function plotSignatures will infer the type of mutations
#and use the correct plot function
plotSignatures(pipeline_subs_res$catalogues,
               output_file = "SNV_catalogues.pdf",
               ncolumns = 2)

#save the annotated mutations
writeTable(pipeline_subs_res$annotated_mutations,"annotated_SNVs.tsv")

#the function plotFitResults can plot both Fit and FitMS objects
plotFitResults(pipeline_subs_res$fitResults,outdir = "signatureFit/")


