#set directory to this file location

#import the package
library(signature.tools.lib)

# set an output directory
outdir <- "~/Example05/"
dir.create(outdir,showWarnings = F,recursive = T)

# load the mutational catalogues
# these are simulated data with 100 catalogues (9 common and 2 rare signatures)
catalogues <- readTable("../../tests/testthat/rareExtraction/SDExample05/catalogues.tsv")

# we begin by exploring the data using clustering
resCl <- cataloguesClustering(catalogues,nclusters = 1:15,
                              outdir = paste0(outdir,"cataloguesClustering/"))

# Let's try to remove samples that may have a rare signatures first
nclustersSelection <- "5"
clusters_to_keep <- c(1:4)
samplesSelected <- rownames(resCl$clusters_table)[resCl$clusters_table[,nclustersSelection] %in% clusters_to_keep]
#now extract common signatures
SignatureExtraction(cat = catalogues[,samplesSelected],
                    outFilePath = paste0(outdir,"Extraction/"),
                    nrepeats = 25,nboots = 4,filterBestOfEachBootstrap = T,
                    nparallel = 4,nsig = 8:11,plotResultsFromAllClusteringMethods = F,
                    parallel = T)

# load the optimal signatures (9)
estimated_signatures <- readTable(paste0(outdir,"Extraction/round_1/sig_9/Sigs_plot_extraction_ns9_nboots4.tsv"))

# is there any sample that was not fully explained by the signatures we found?
unexplSamples <- unexplainedSamples(outfileRoot = paste0(outdir,"unexplained/Example05"),
                                    catalogues = catalogues,
                                    sigs = estimated_signatures,
                                    nmuts_threshold = 300,
                                    pvalue_threshold = 0.15)
# yes, 4 samples, get their residual
significant_residuals <- unexplSamples$all_residuals[,unexplSamples$which_significant]
# now we cluster the residuals
resCl_residuals <- cataloguesClustering(significant_residuals,
                                        nclusters = 1:3,
                                        outdir = paste0(outdir,"residualClustering/"))

# we decide that there are two clusters, having two samples in each
# now we need to extract the rare signatures
resRareSigs <- rareSignatureExtraction(outfileRoot = paste0(outdir,"ExtractionRare/Example05"),
                                       catalogues = catalogues,
                                       residuals = unexplSamples$all_residuals,
                                       unexpl_samples = unexplSamples$unexplSamples,
                                       clusters = resCl_residuals$clusters_table[,"2"],
                                       useclusters = list(c(1),c(2)),
                                       commonSignatures = estimated_signatures,
                                       commonExposures = unexplSamples$exposures)

# now let's get the exposures
resFinalExpo <- finaliseCommonRareSignatureExposures(outfileRoot = paste0(outdir,"ExtractionRare/Example05"),
                                                     catalogues = catalogues,
                                                     commonSigs = estimated_signatures,
                                                     listofsignatures = resRareSigs$listofsignatures,
                                                     listofsamples = resRareSigs$listofsamples,
                                                     nboot = 50,nparallel = 4)



