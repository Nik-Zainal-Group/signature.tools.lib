#!/usr/bin/env Rscript

library(getopt)

how_to <- function(){
  message(" ")
  message("This script runs the HRDetect pipeline.")
  message("Run this script as follows:")
  message(" ")
  message("hrDetect [-i INPUTTABLE] [-o OUTDIR] [-b] [-s SIGTYPE] [-O organ] [-e GENOMEV] [-n NPARALLEL] [-f NBOOTFIT] [-r RANDOMSEED]")
  message(" ")
  message("    -i INPUTTABLE    Tab separate input table with the list of files for each sample. Columns should be: sample, SNV_vcf_files, SNV_tab_files, Indels_vcf_files, Indels_tab_files, CNV_tab_files, SV_bedpe_files. Note that only one column of SNV_vcf_files and SNV_tab_files is necessary")
  message("    -o OUTDIR        Name of the output directory. If omitted a name will be given automatically.")
  message("    -b               Request HRDetect with bootstrap")
  message("    -s SIGTYPE       Either COSMICv2, COSMICv3.2, RefSigv1 or RefSigv2. When selecting RefSigv2, signature fitting for SNVs will be performed with FitMS")
  message("    -O ORGAN         When using RefSigv1 or RefSigv2 as SIGTYPE you need to specify the organ or your samples, as organ-specific signatures will be used. Use one of the following organs: Biliary, Bladder, Bone_SoftTissue, Breast, Cervix (v1 only), CNS, Colorectal, Esophagus, Head_neck, Kidney, Liver, Lung, Lymphoid, NET (v2 only), Oral_Oropharyngeal (v2 only), Ovary, Pancreas, Prostate, Skin, Stomach, Uterus")
  message("    -l COSMICLIST    If SIGTYPE is COSMICv2, specify a comma separated list of number of cosmic signatures to use, such as: 1,2,5,8,13. If SIGTYPE is COSMICv3.2, specify the SBS names instead, such as: SBS1,SBS7a,SBS10a,SBS18. If no list is specified, all signatures will be used")
  message("    -e GENOMEV       Indicate the genome version to be used: hg19, hg38 or mm10")
  message("    -n NPARALLEL     Number of parallel CPUs to be used")
  message("    -f NBOOTFIT      Number of bootstrap to be used in signature fit, default is 100")
  message("    -r RANDOMSEED    Specify a random seed to obtain always the same identical results")
  message("    -h               Show this explanation")
  message(" ")
}

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'input',         'i', 1, "character",
  'help',          'h', 0, "logical",
  'bootstrap',     'b', 0, "logical",
  'outdir',        'o', 1, "character",
  'sigtype',       's', 1, "character",
  'organ',         'O', 1, "character",
  'cosmiclist',    'l', 1, "character",
  'nparallel',     'n', 1, "double",
  'genomev',       'e', 1, "character",
  'nbootFit',      'f', 1, "double",
  'randomSeed',    'r', 1, "double"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# if help was asked for print a friendly message
# and exit with a zero error code
if ( !is.null(opt$help) ) {
  how_to()
  q(status=0,save = "no")
}

# set some reasonable defaults for the options that are needed,
# but were not specified.
bootstrap_scores <- FALSE
cosmic_siglist <- NULL
organ <- NULL

if ( !is.null(opt$bootstrap) ) {
  bootstrap_scores <- TRUE
}

if ( is.null(opt$input ) ) { 
  message("\nMissing INPUTTABLE. Quit.\n")
  how_to()
  q(status=1,save = "no")
}
if ( is.null(opt$outdir) ) { 
  opt$outdir = "hrDetectResults/"    
}
if ( is.null(opt$sigtype) ) { 
  opt$sigtype = "COSMICv2"    
}
if ( !is.null(opt$cosmiclist) ) {
  if(opt$sigtype == "COSMICv2"){
    cosmic_siglist <- as.numeric(unlist(strsplit(opt$cosmiclist,split = ",")))
  }else if (opt$sigtype == "COSMICv3.2"){
    cosmic_siglist <- unlist(strsplit(opt$cosmiclist,split = ","))
  }
}
if ( is.null(opt$nparallel) ) { 
  opt$nparallel = 1   
}
if ( is.null(opt$nbootFit) ) { 
  opt$nbootFit = 100  
}
if ( is.null(opt$genomev ) ) { 
  message("\nMissing GENOMEV. Quit.\n")
  how_to()
  q(status=1,save = "no")
}

library(signature.tools.lib)

input_table <- opt$input
outdir <- opt$outdir
genomev <- opt$genomev
signature_type <- opt$sigtype
nparallel <- opt$nparallel
organ <- opt$organ
nboots <- opt$nbootFit
randomSeed <- opt$randomSeed

dir.create(outdir,showWarnings = FALSE,recursive = TRUE)

#retrieve samples and file names to be used from the files table
input_t <- read.table(input_table,sep = "\t",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)

#initialise input variables
samples <- NULL
SNV_vcf_files <- NULL
SNV_tab_files <- NULL
Indels_vcf_files <- NULL
Indels_tab_files <- NULL
CNV_tab_files <- NULL
SV_bedpe_files <- NULL

#read input file lists
if ("sample" %in% colnames(input_t)) {
  samples <- input_t$sample
}
if ("SNV_vcf_files" %in% colnames(input_t)) {
  SNV_vcf_files <- input_t$SNV_vcf_files
  names(SNV_vcf_files) <- samples
}
if ("SNV_tab_files" %in% colnames(input_t)) {
  SNV_tab_files <- input_t$SNV_tab_files
  names(SNV_tab_files) <- samples
}
if ("Indels_vcf_files" %in% colnames(input_t)) {
  Indels_vcf_files <- input_t$Indels_vcf_files
  names(Indels_vcf_files) <- samples
}
if ("Indels_tab_files" %in% colnames(input_t)) {
  Indels_tab_files <- input_t$Indels_tab_files
  names(Indels_tab_files) <- samples
}
if ("CNV_tab_files" %in% colnames(input_t)) {
  CNV_tab_files <- input_t$CNV_tab_files
  names(CNV_tab_files) <- samples
}
if ("SV_bedpe_files" %in% colnames(input_t)) {
  SV_bedpe_files <- input_t$SV_bedpe_files
  names(SV_bedpe_files) <- samples
}

#run HRDetect
hrdet_res <- HRDetect_pipeline(genome.v = genomev,
                               SNV_vcf_files = SNV_vcf_files,
                               SNV_tab_files = SNV_tab_files,
                               Indels_vcf_files = Indels_vcf_files,
                               Indels_tab_files = Indels_tab_files,
                               CNV_tab_files = CNV_tab_files,
                               SV_bedpe_files = SV_bedpe_files, 
                               bootstrapHRDetectScores = bootstrap_scores,
                               signature_type = signature_type,
                               organ = organ,
                               cosmic_siglist=cosmic_siglist,
                               nbootFit = nboots,
                               randomSeed = randomSeed,
                               nparallel = nparallel)

#save object
save(hrdet_res,file = paste0(outdir,"/hrdet_res.rData"))

#saving tables
message("Writing table data_matrix...")
write.table(hrdet_res$data_matrix,sep = "\t",
            file = paste0(outdir,"/data_matrix.tsv"),
            col.names = TRUE,row.names = TRUE,quote = FALSE)
message("Writing table hrdetect_output...")
write.table(hrdet_res$hrdetect_output,sep = "\t",
            file = paste0(outdir,"/hrdetect_output.tsv"),
            col.names = TRUE,row.names = TRUE,quote = FALSE)
if (bootstrap_scores) {
  message("Writing table hrdetect_output_confidence...")
  write.table(hrdet_res$q_5_50_95,sep = "\t",
                                  file = paste0(outdir,"/hrdetect_output_confidence.tsv"),
                                  col.names = TRUE,row.names = TRUE,quote = FALSE)
}

message("Writing SNV and SV catalogues tables...")
write.table(hrdet_res$SNV_catalogues,sep = "\t",
            file = paste0(outdir,"/SNV_catalogues.tsv"),
            col.names = TRUE,row.names = TRUE,quote = FALSE)
write.table(hrdet_res$SV_catalogues,sep = "\t",
            file = paste0(outdir,"/SV_catalogues.tsv"),
            col.names = TRUE,row.names = TRUE,quote = FALSE)

message("Plotting SNV and SV catalogues...")
plotSignatures(signature_data_matrix = hrdet_res$SNV_catalogues,
               output_file = paste0(outdir,"/SNV_catalogues.pdf"))
plotSignatures(signature_data_matrix = hrdet_res$SV_catalogues,
               output_file = paste0(outdir,"/SV_catalogues.pdf"))

message("Writing siganture exposure tables...")
write.table(hrdet_res$exposures_subs,sep = "\t",
            file = paste0(outdir,"/Exposures_subs.tsv"),
            col.names = TRUE,row.names = TRUE,quote = FALSE)
write.table(hrdet_res$exposures_rearr,sep = "\t",
            file = paste0(outdir,"/Exposures_rearr.tsv"),
            col.names = TRUE,row.names = TRUE,quote = FALSE)

message("Writing table indels_classification_table...")
write.table(hrdet_res$indels_classification_table,sep = "\t",
            file = paste0(outdir,"/indels_classification_table.tsv"),
            col.names = TRUE,row.names = FALSE,quote = FALSE)

#plotting
message("Plotting HRDetect contributions...")
plot_HRDetect_Contributions(file_name = paste0(outdir,"/HRDetect_contributions.pdf"),
                                   hrdetect_output = hrdet_res$hrdetect_output)

message("Plotting HRDetect overall plot...")
plot_HRDetect_overall(file_name = paste0(outdir,"/HRDetect_overall.jpg"),
                      hrdetect_output = hrdet_res$hrdetect_output)

if (bootstrap_scores) {
  message("Plotting HRDetect with Bootstrap...")
  plot_HRDetect_BootstrapScores(outdir = outdir,
                                hrdetect_res = hrdet_res)
}

message("Plotting Subs and Rearr Signature Fit files...")
if(signature_type=="RefSigv2"){
  plotFitMS(hrdet_res$fitRes_subs,outdir = paste0(outdir,"/subsfit/"))
}else{
  plotFit(hrdet_res$fitRes_subs,outdir = paste0(outdir,"/subsfit/"))
}
plotFit(hrdet_res$fitRes_rearr,outdir = paste0(outdir,"/rearrfit/"))

#genomeplots
dir.create(paste0(outdir,"/genomeplots/"),showWarnings = FALSE,recursive = TRUE)
for (s in samples) {
  SNV_vcf_file <- SNV_vcf_files[s]
  if(!is.null(SNV_vcf_file)){
    if(is.na(SNV_vcf_files[s])) SNV_vcf_file <- NULL
  }
  if(!is.null(SNV_vcf_file)){
    message("Generating genomeplot for ",s,"...")
    genomePlot(subsVcf.file = SNV_vcf_file,
               indelsVcf.file = Indels_vcf_files[s],
               cnvsTab.file = CNV_tab_files[s],
               rearrBedpe.file = SV_bedpe_files[s],
               sampleID = s,
               genome.v = genomev,
               out_path = paste0(outdir,"/genomeplots/"))
  }else{
    message("Skipping genomeplot for ",s," because of missing SNV VCF file.")
  }
}
message(" ")
message("Done!")
message(" ")

q(status = 0,save = "no")
