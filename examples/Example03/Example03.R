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

#load SNV data and convert to SNV mutational catalogues
SNVcat_list <- list()
for (i in 1:length(SNV_tab_files)){
  tmpSNVtab <- read.table(SNV_tab_files[i],sep = "\t",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
  #convert to SNV catalogue, see ?tabToSNVcatalogue or ?vcfToSNVcatalogue for details
  res <- tabToSNVcatalogue(subs = tmpSNVtab,genome.v = "hg19")
  colnames(res$catalogue) <- sample_names[i]
  SNVcat_list[[i]] <- res$catalogue
}
#bind the catalogues in one table
SNV_catalogues <- do.call(cbind,SNVcat_list)

#the catalogues can be plotted as follows
plotSubsSignatures(signature_data_matrix = SNV_catalogues,plot_sum = TRUE,output_file = "SNV_catalogues.pdf")

#perform signature fit using a multi-step approach where organ-specific common and rare signatures are used
subs_fit_res <- FitMS(catalogues = SNV_catalogues,
                      exposureFilterType = "giniScaledThreshold",
                      useBootstrap = TRUE,
                      organ = "Breast")
plotFitMS(subs_fit_res,outdir = "signatureFit/")

#The signature exposures can be found here and correspond to the median of the boostrapped runs followed by false positive filters. See ?Fit for details
snv_exp <- subs_fit_res$exposures

#write the results
writeTable(snv_exp,"RefSigSubsExposures.tsv")

