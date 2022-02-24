#!/usr/bin/env Rscript

library(getopt)

how_to <- function(){
  message(" ")
  message("This script runs the signature fit pipeline of the R package")
  message("signatures.tools.lib, using the Fit or FitMS functions.")
  message(" ")
  message("Run this script as follows:")
  message(" ")
  message("signatureFit [OPTIONS]")
  message(" ")
  message("Available options:")
  message("  -o, --outdir=DIR          Name of the output directory. If omitted a name")
  message("                              will be given automatically.")
  message("  -b, --bootstrap           Request signature fit with bootstrap")
  message("  -x, --snvvcf=SNVVCF       SNVVCF is a tab separated file containing two")
  message("                              columns. The first column contains the sample")
  message("                              names, while the second column contains the")
  message("                              corresponding SNV vcf file names.")
  message("  -X, --snvtab=SNVTAB       SNVTAB is a tab separated file containing two")
  message("                              columns. The first column contains the sample")
  message("                              names, while the second column contains the")
  message("                              corresponding SNV tab file names. Each SNV tab")
  message("                              file should have a header with the following")
  message("                              columns: chr, position, REF, ALT.")
  message("  -s, --sigversion=SIGVERSION")
  message("                            Either COSMICv2, COSMICv3.2, RefSigv1 or RefSigv2.")
  message("                              If not specified SIGVERSION=RefSigv2.")
  message("  -O, --organ=ORGAN         When using RefSigv1 or RefSigv2 as SIGVERSION,")
  message("                              organ-specific signatures will be used.")
  message("                              If SIGVERSION is COSMICv2 or COSMICv3.2, then a")
  message("                              selection of signatures found in the given organ")
  message("                              will be used. Organ names depend on the selected")
  message("                              SIGVERSION. For RefSigv1 or RefSigv2: Biliary,")
  message("                              Bladder, Bone_SoftTissue, Breast, Cervix (v1 only),")
  message("                              CNS, Colorectal, Esophagus, Head_neck, Kidney,")
  message("                              Liver, Lung, Lymphoid, NET (v2 only),")
  message("                              Oral_Oropharyngeal (v2 only), Ovary, Pancreas,")
  message("                              Prostate, Skin, Stomach, Uterus.")
  message("  -l, --signames=SIGNAMES   If no ORGAN is specified, SIGNAMES can be used to")
  message("                              provide a comma separated list of signature names")
  message("                              to select from the COSMIC or reference signatures,")
  message("                              depending on the SIGVERSION requested. For example,")
  message("                              for COSMICv3.2 use: SBS1,SBS2,SBS3.")
  message("  -e, --genomev=GENOMEV     Genome version to be used: hg19, hg38 or mm10.")
  message("                              If not specified GENOMEV=hg19.")
  message("  -m, --fitmethod=FITMETHOD Either Fit or FitMS. If not specified FITMETHOD=FitMS")
  message("  -n, --nparallel=NPARALLEL Number of parallel CPUs to be used.")
  message("  -f, --nboot=NBOOT         Number of bootstrap to be used when bootstrap is")
  message("                              requested (-b), if not specified, NBOOT=200.")
  message("  -r, --randomSeed=SEED     Specify a random seed to obtain always the same")
  message("                              identical results.")
  message("  -h, --help                Show this explanation.")
  message(" ")
}


# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'help',          'h', 0, "logical",
  'bootstrap',     'b', 0, "logical",
  'outdir',        'o', 1, "character",
  'snvvcf',        'x', 1, "character",
  'snvtab',        'X', 1, "character",
  'sigversion',    's', 1, "character",
  'organ',         'O', 1, "character",
  'signames',      'l', 1, "character",
  'fitmethod',     'm', 1, "character",
  'nparallel',     'n', 1, "double",
  'genomev',       'e', 1, "character",
  'nboot',         'f', 1, "double",
  'randomSeed',    'r', 1, "double"
), byrow=TRUE, ncol=4)

# collect the options
opt = getopt(spec)

# check options

# if help was asked for, print a friendly message and quit
if ( !is.null(opt$help) ) {
  how_to()
  q(status=0,save = "no")
}

if ( !is.null(opt$bootstrap) ) {
  useBootstrap <- TRUE
}else{
  useBootstrap <- FALSE
}

if ( !is.null(opt$outdir) ) {
  outdir <- opt$outdir
  if(substr(outdir,nchar(outdir),nchar(outdir))!="/" & substr(outdir,nchar(outdir),nchar(outdir))!="\\"){
    outdir <- paste0(outdir,"/")
  }
}else{
  outdir <- "signatureFitResults/"
}

# Load location of SNV vcf files
if(!is.null(opt$snvvcf)){
  if(file.exists(opt$snvvcf)){
    infotable <- read.table(opt$snvvcf,sep = "\t",stringsAsFactors = FALSE)
    snvvcf <- infotable[,2]
    names(snvvcf) <- infotable[,1]
  }else{
    message("snvvcf file does not exist: ",opt$snvvcf)
    q(status=1,save = "no")
  }
}else{
  snvvcf <- NULL
}

# Load location of SNV tab files
if(!is.null(opt$snvtab)){
  if(file.exists(opt$snvtab)){
    infotable <- read.table(opt$snvtab,sep = "\t",stringsAsFactors = FALSE)
    snvtab <- infotable[,2]
    names(snvtab) <- infotable[,1]
  }else{
    message("snvtab file does not exist: ",opt$snvtab)
    q(status=1,save = "no")
  }
}else{
  snvtab <- NULL
}

if ( !is.null(opt$sigversion) ) {
  sigversion <- opt$sigversion
}else{
  sigversion <- "RefSigv2"
}

if ( !is.null(opt$signames) ) {
  signames <- unlist(strsplit(opt$signames,split = ","))
}else{
  signames <- NULL
}

if ( !is.null(opt$fitmethod) ) {
  fitmethod <- opt$fitmethod
}else{
  fitmethod <- "FitMS"
}

if ( !is.null(opt$genomev) ) {
  genomev <- opt$genomev
}else{
  genomev <- "hg19"
}

if ( !is.null(opt$nparallel) ) {
  nparallel <- opt$nparallel
}else{
  nparallel <- 1
}

if ( !is.null(opt$nboot) ) {
  nboots <- opt$nboot
}else{
  nboots <- 200
}

organ <- opt$organ
randomSeed <- opt$randomSeed

# create output directory
dir.create(outdir,showWarnings = F,recursive = T)

# run the pipeline
pipelineRes <- signature.tools.lib::signatureFit_pipeline(catalogues=NULL,
                                  genome.v=genomev,
                                  organ=organ,
                                  SNV_vcf_files=snvvcf,
                                  SNV_tab_files=snvtab,
                                  DNV_vcf_files=NULL,
                                  DNV_tab_files=NULL,
                                  SV_bedpe_files=NULL,
                                  signatures=NULL,
                                  rare_signatures=NULL,
                                  signature_version=sigversion,
                                  signature_names=signames,
                                  fit_method=fitmethod, 
                                  optimisation_method = "KLD",
                                  useBootstrap = useBootstrap,
                                  nboot = nboots,
                                  exposureFilterType = "fixedThreshold", # or "giniScaledThreshold"
                                  threshold_percent = 5,
                                  giniThresholdScaling = 10,
                                  multiStepMode = "errorReduction", # or "partialNMF", or "errorReduction", or "cossimIncrease"
                                  threshold_p.value = 0.05,
                                  rareSignatureTier = "T2",  #either T1 for observed in organ or T2 for extended
                                  residualNegativeProp = 0.003,
                                  minResidualMutations = NULL,
                                  minCosSimRareSig = 0.8,
                                  minErrorReductionPerc = 15,
                                  minCosSimIncrease = 0.02,
                                  maxRareSigsPerSample = 1,
                                  nparallel = nparallel,
                                  randomSeed = randomSeed,
                                  verbose = FALSE)

if(!is.null(pipelineRes$fitResults)){
  # save fit to file
  signature.tools.lib::saveFitToFile(fitObj = pipelineRes$fitResults,
                                     filename = paste0(outdir,"fitData.rData"))
  # plot fit results
  signature.tools.lib::plotFitResults(fitObj = pipelineRes$fitResults,
                                      outdir = outdir)
}


q(status = 0,save = "no")
