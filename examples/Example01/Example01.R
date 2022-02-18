#set directory to this file location

#import the package
library(signature.tools.lib)

#set sample names
sample_names <- c("sample1","sample2")

#set the file names. 
SNV_tab_files <- c("../../tests/testthat/test_hrdetect_1/test_hrdetect_1.snv.simple.txt",
                   "../../tests/testthat/test_hrdetect_2/test_hrdetect_2.snv.simple.txt")
SV_bedpe_files <- c("../../tests/testthat/test_hrdetect_1/test_hrdetect_1.sv.bedpe",
                    "../../tests/testthat/test_hrdetect_2/test_hrdetect_2.sv.bedpe")
Indels_vcf_files <- c("../../tests/testthat/test_hrdetect_1/test_hrdetect_1.indel.vcf.gz",
                      "../../tests/testthat/test_hrdetect_2/test_hrdetect_2.indel.vcf.gz")
CNV_tab_files <- c("../../tests/testthat/test_hrdetect_1/test_hrdetect_1.cna.txt",
                   "../../tests/testthat/test_hrdetect_2/test_hrdetect_2.cna.txt")

#name the vectors entries with the sample names
names(SNV_tab_files) <- sample_names
names(SV_bedpe_files) <- sample_names
names(Indels_vcf_files) <- sample_names
names(CNV_tab_files) <- sample_names

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

#fit the 12 breast cancer signatures using the bootstrap signature fit approach
sigsToUse <- c(1,2,3,5,6,8,13,17,18,20,26,30)
subs_fit_res <- Fit(catalogues = SNV_catalogues,
                    signatures = COSMIC30_subs_signatures[,sigsToUse],
                    useBootstrap = TRUE,
                    nboot = 100,
                    nparallel = 4)
plotFit(subs_fit_res,outdir = "signatureFit/")

#The signature exposures can be found here and correspond to the median of the boostrapped runs followed by false positive filters. See ?Fit for details
snv_exp <- subs_fit_res$exposures

#The HRDetect pipeline will compute the HRDetect probability score for the samples to be Homologous Recombination Deficient. HRDetect is a logistic regression classifier that requires 6 features to compute the probability score. These features can be supplied directly in an input matrix, or pipeline can compute these features for you if you supply the file names. It is possible to supply a mix of features and file names.

#Initialise feature matrix
col_hrdetect <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
input_matrix <- matrix(NA,nrow = length(sample_names),ncol = length(col_hrdetect),dimnames = list(sample_names,col_hrdetect))

#We have already quantified the amount of SNV signatures in the samples, so we can supply these via the input matrix
input_matrix[rownames(snv_exp),"SNV3"] <- snv_exp[,"Signature3"]
input_matrix[rownames(snv_exp),"SNV8"] <- snv_exp[,"Signature8"]

#run the HRDetect pipeline, for more information see ?HRDetect_pipeline
res <- HRDetect_pipeline(input_matrix,
                         genome.v = "hg19",
                         signature_type = "COSMICv2",
                         SV_bedpe_files = SV_bedpe_files,
                         Indels_vcf_files = Indels_vcf_files,
                         CNV_tab_files = CNV_tab_files,
                         nparallel = 2)

#save HRDetect scores
writeTable(res$hrdetect_output,file = "HRDetect_res.tsv")

