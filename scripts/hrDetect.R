#!/usr/bin/env Rscript

library(getopt)

how_to <- function(){
  message(" ")
  message("This script runs the HRDetect pipeline.")
  message("Run this script as follows:")
  message(" ")
  message("hrDetect [-i INPUTTABLE] [-o OUTDIR] [-b] [-s SIGTYPE] [-e GENOMEV] [-n NPARALLEL]")
  message(" ")
  message("    -i INPUTTABLE    Tab separate input table with the list of files for each sample. Columns should be: sample, SNV_vcf_files, SNV_tab_files, Indels_vcf_files, CNV_tab_files, SV_bedpe_files. Note that only one column of SNV_vcf_files and SNV_tab_files is necessary")
  message("    -o OUTDIR        Name of the output directory. If omitted a name will be given automatically.")
  message("    -b               Request HRDetect with bootstrap")
  message("    -s SIGTYPE       Either COSMIC or one of the following organs: Biliary, Bladder, Bone_SoftTissue, Breast, Cervix, CNS, Colorectal, Esophagus, Head_neck, Kidney, Liver, Lung, Lymphoid, Ovary, Pancreas, Prostate, Skin, Stomach, Uterus")
  message("    -l COSMICLIST    If SIGTYPE is COSMIC, specify a comma separated list of number of cosmic signatures to use, such as: 1,2,5,8,13")
  message("    -e GENOMEV       Indicate the genome version to be used: hg19, hg38 or mm10")
  message("    -n NPARALLEL     Number of parallel CPUs to be used")
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
  'cosmiclist',    'l', 1, "character",
  'nparallel',     'n', 1, "double",
  'genomev',       'e', 1, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  how_to()
  q(status=1,save = "no")
}

# set some reasonable defaults for the options that are needed,
# but were not specified.
bootstrap_scores <- FALSE

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
  opt$sigtype = "COSMIC"    
}
if ( !is.null(opt$cosmiclist) ) {
  cosmic_siglist <- as.numeric(unlist(strsplit(opt$cosmiclist,split = ",")))
}
if ( is.null(opt$nparallel) ) { 
  opt$nparallel = 1   
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

dir.create(outdir,showWarnings = FALSE,recursive = TRUE)

#retrieve samples and file names to be used from the files table
input_t <- read.table(input_table,sep = "\t",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)

#initialise input variables
samples <- NULL
SNV_vcf_files <- NULL
SNV_tab_files <- NULL
Indels_vcf_files <- NULL
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
                               CNV_tab_files = CNV_tab_files,
                               SV_bedpe_files = SV_bedpe_files,
                               bootstrap_scores = bootstrap_scores,
                               signature_type = signature_type,
                               cosmic_siglist=cosmic_siglist,
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
plot_HRDetect_Contributions(file_name = paste0(outdir,"/HRDetect_contributions.jpg"),
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
plot_SignatureFit_withBootstrap(outdir = paste0(outdir,"/subsfit/"),
                                boostrapFit_res = hrdet_res$bootstrap_fit_subs,
                                type_of_mutations = "subs")
plot_SignatureFit_withBootstrap(outdir = paste0(outdir,"/rearrfit/"),
                                boostrapFit_res = hrdet_res$bootstrap_fit_rearr,
                                type_of_mutations = "rearr")

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
