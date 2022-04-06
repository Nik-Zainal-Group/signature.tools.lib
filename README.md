# signature.tools.lib R package

## Table of content

- [Introduction to the package](#intro)
- [Versions](#version)
- [How to cite us](#cite)
- [Systems Requirements](#req)
- [Testing the package](#test)
- [How to use this package](#howtouse)
- [Package documentation](#docs)
- [Functions provided by the package](#functions)
- [Command line scripts](#scripts)
- [Examples](#examples)
  - [Test examples](#examplestests)
  - [Example 01](#examplese01)
  - [Example 02](#examplese02)
  - [Example 03](#examplese03)
  - [Example 04](#examplese04)
- [Frequently Asked Questions](#faq)

<a name="intro"/>

## Introduction to the package

```signature.tools.lib``` is an R package for mutational signatures analysis.
Mutational signatures are patterns of somatic mutations that reveal what
mutational processes have been active in a cell. Mutational processes
can be caused by exposure to mutagens, such as chemicals present in
cigarettes, or defects in DNA repair pathways, such as homologous
recombination repair.

The package supports hg19 and hg38 as well as mm10. It provides our
latest algorithms for signature fit and extraction, as well as various
utility functions and the HRDetect pipeline. The list and description of
the most important functions is given below.

## Versions

<a name="version"/>

2.1.2

- Signature fit objects obtained from Fit or FitMS can be passed to the HRDetect pipeline

2.1.1

- Updated the ```HRDetect_pipeline``` function to use the ```signatureFit_pipeline``` function

2.1.0

- Added the ```signatureFit_pipeline``` function, which is a flexible interface for signature fit analysis
- Added the ```signatureFit``` command line script, which is a wrapper for the ```signatureFit_pipeline``` function

2.0.1

- Updated HRDetect pipeline to work with FitMS
- Updated DNV catalogues functions

2.0

- New signature fit multi-step algorithm, FitMS, from Degasperi et al. 2022, *Science*
- New organ-specific signatures obtained using Genomics England Cancer data, from Degasperi et al. 2022, *Science*
- Reorganisation of old signature fit functions under the new Fit function
- Rewrite of plot scripts for displaying signature fit results, now available using plotFit and plotFitMS functions
- Added support for trinucleotide variant catalogues
- New export functions for converting Fit and FitMS results to JSON

1.0

- Signature analysis functions for signature extraction and fit
- Mutational signatures available are COSMICv2 and reference signatures from Degasperi et al. 2020, *Nature Cancer*
- HRDetect pipeline available as command line script, includes HRDetect with bootstrap as published in Degasperi et al. 2020, *Nature Cancer*
- Genome plot circos function (fork of Wellcome Sanger Institute genome plot function)
- HRD indexes functions written by Nicolai Juul Birkbak

## How to cite us

<a name="cite"/>

These are the research articles associated with ```signature.tools.lib```:

- A. Degasperi et al. **Substitution mutational signatures in whole-genome-sequenced cancers in the UK population.** *Science*, doi:10.1126/science.abl9283, 2022.
- A. Degasperi et al. **A practical framework and online tool for mutational signature analyses show intertissue variation and driver dependencies.** *Nature Cancer*, doi:10.1038/s43018-020-0027-5, 2020.

More in details, these are the specific signatures and algorithms introduced in each publication:

- Degasperi et al. 2022, *Science*: RefSig SBS v2, RefSig DBS v1, FitMS
- Degasperi et al. 2020, *Nature Cancer*: RefSig SBS v1, RefSig Rearrangements (SV) v1, Fit with bootstrap, HRDetect with bootstrap.

<a name="req"/>

## Systems Requirements

No special hardware is required to run this software. This is an R
package so it will work on any computer with R installed. We recommend to
use the latest version of R, as some of the R dependencies constantly increase
the minimum version of R required for installation.

You can install ```signature.tools.lib``` by entering the R environment from the main
directory and typing:

```
install.packages("devtools")
devtools::install()
```

The above commands will take care of all the dependencies. In some
cases, it may happen that some dependencies cannot be installed
automatically. In this case, an error will indicate which package could
not be found, then simply install these packages manually. The
installation should not take more than a few minutes.

This is the full list of R package dependencies:

```
    VariantAnnotation,
    BSgenome.Hsapiens.UCSC.hg38,
    BSgenome.Hsapiens.1000genomes.hs37d5,
    BSgenome.Mmusculus.UCSC.mm10,
    BSgenome.Cfamiliaris.UCSC.canFam3,
    SummarizedExperiment,
    BiocGenerics,
    GenomeInfoDb,
    NMF,
    foreach,
    doParallel,
    lpSolve,
    ggplot2,
    methods,
    cluster,
    stats,
    NNLM,
    nnls,
    GenSA,
    gmp,
    plyr,
    RCircos,
    scales,
    GenomicRanges,
    IRanges,
    BSgenome, 
    readr,
    doRNG,
    combinat
```

We have noticed that the ```NNLM``` package is frequently unavailable to download automatically
during the R installation, so we have changed the installation process to install ```NNLM``` from
the github repository [linxihui/NNLM](https://github.com/linxihui/NNLM).

<a name="test"/>

## Testing the package

You can test the package by entering the package main directory and
typing from the R environment:

```
devtools::test()
```

<a name="howtouse"/>

## How to use this package

**PLEASE NOTE:** project-specific file conversions and filtering should
be done before and outside the use of this library. This library should
only report to the user the data format errors and suggest how to
correct them, while it is the responsibility of the user to supply the
necessary data formatted correctly with the necessary features and
correct column names.

<a name="docs"/>

## Package documentation

**DOCUMENTATION:** Documentation for each of the functions below is
provided as R documentation, and it is installed along with the R
package. The documentation should give detailed explanation of the input
data required, such as a list of data frame columns and their
explanation. To access the documentation you can use the ```?function```
syntax in R, for each of the functions below. For example, in R or
RStudio, type ```?HRDetect_pipeline```.

<a name="functions"/>

## Functions provided by the package

Functions for single nucleotide variants:

- **```vcfToSNVcatalogue(...)```**: converts a VCF file into a SNV
catalogue.
- **```tabToSNVcatalogue(...)```**: converts a data frame into a SNV
catalogue. The data frame should have the following minimal columns:
chr, position, REF, ALT.
- **```plotSubsSignatures(...)```**: function to plot one or more
substitution signatures or catalogues.

Functions for rearrangements:

- **```bedpeToRearrCatalogue(...)```**: converts a data frame (loaded
from a BEDPE file) into a rearrangement catalogue. If the columns
```is.clustered``` or ```svclass``` are not present, this function will
attempt to compute them.
- **```plotRearrSignatures(...)```**: function to plot one or more
rearrangement signatures or catalogues.

Functions for dinucleotide variants:

- **```vcfToDNVcatalogue(...)```**: given a vcf file containing a list of single and double nucleotide variants (SNVs and DNVs),
this function initially attempts to find all possible DNVs by checking adjacent SNVs, and then
annotates the DNVs and finally generates a DNV catalogue. The catalogue channels are in Zou's style.
- **```tabToDNVcatalogue(...)```**: same as ```vcfToDNVcatalogue(...)``` but takes a data frame as input.
- **```convertToAlexandrovChannels(...)```**: Function to convert DNV signatures or catalogues
from Zou's to Alexandrov's style. 
- **```plotDNVSignatures(...)```**: plot one or more DNV signatures or catalogues, compatible with both
Zou's style and Alexandrov's style of channels.

Functions for signature extraction and signature fit:

- **```SignatureExtraction(...)```**: perform signature extraction, by
applying NMF to the input matrix. Multiple NMF runs and bootstrapping is
used for robustness, followed by clustering of the solutions. A range of
number of signatures to be used is required.
- **```Fit(...)```**: This is a standard interface for basic signature fit with/without bootstrap.
The object returned by this function can be passed to the ```plotFit()``` function for automated plotting of the results. Use the function
```fitToJSON``` to export the results into a JSON file.
- **```FitMS(...)```**: Given a set of mutational catalogues, this function will attempt fit mutational signature in a multi-step manner. 
In the first step, only a set of common signatures are fitted into the samples. In the following steps, one or more rare signatures
are fitted into the samples in addition to the common signatures. Common and rare signatures can be determined automatically 
by providing the name of an organ, or can be supplied by the user.
The object returned by this function can be passed to the ```plotFitMS()``` function for automated plotting of the results. Use the function
```fitMStoJSON``` to export the results into a JSON file.
A manual for FitMS can be found in the ```userManuals``` folder.
- **```signatureFit_pipeline```**: an interface for the ```Fit``` and ```FitMS``` functions, which aim to automate various signature fit analysis steps, like generating the mutational catalogues and selecting which mutational signatures to fit. This function can be accessed via command line using the ```signatureFit``` script in the ```scripts``` folder.
- **```plotFitResults(...)```**: this function can be used to plot results objects from both ```Fit``` and ```FitMS``` functions. The object type will be inferred automatically and either ```plotFit()``` or ```plotFitMS()``` will be used.

Functions for organ-specific signatures and exposures conversion:

- **```getOrganSignatures(...)```**: returns organ-specific signatures of a requested organ and mutation type.
- **```convertExposuresFromOrganToRefSigs(...)```**: after signature fit using organ-specific signatures,
use this function to apply the conversion matrix, to convert organ-specific signature exposures into reference signature exposures. 

Functions for HRD indexes:

- **```ascatToHRDLOH(...)```**: compute the HRD-LOH index from a data
frame formatted similarly to an ASCAT output file. This is a wrapper for
the function ```calc.hrd```, with a simplified interface, where only a
data frame is requested.
- **```plotCopyNumbers(...)```**: Plot the copy numbers across the Chromosomes.
Optionally, plot also highlight regions at the bottom.
For example, could be used to highlight HRD-LOH regions..
- **```calc.hrd(...)```**: compute HRD-LOH (Loss of Heterozygosity),
written by Nicolai Juul Birkbak
- **```calc.ai(...)```**: compute HRD-NtAI (Number of telomeric Allelic
Imbalances), written by Nicolai Juul Birkbak
- **```calc.lst(...)```**: compute HRD-LST (Large-scale state
transitions), written by Nicolai Juul Birkbak

Functions for indels classification:

- **```vcfToIndelsClassification(...)```**: converts a VCF file
containing indels into a data frame where the indels are classified.
Also returns a summary of count and proportion of the various classes of
indels.
- **```tabToIndelsClassification(...)```**: converts a data frame
containing indels into a data frame where the indels are classified.
Also returns a summary of count and proportion of the various classes of
indels.

Functions for HRDetect:

- **```HRDetect_pipeline(...)```**: this is a flexible interface to the
HRDetect pipeline to compute the HRDetect BRCAness probability score, as
published in Davies et al. 2017.
- **```applyHRDetectDavies2017(...)```**: given a data frame with
samples as rows and features as columns, this function will compute the
HRDetect BRCAness probability for each sample. The following six
features should be included in the matrix: 1) proportion of deletions
with micro-homology, 2) number of mutations of substitution signature 3, 3)
number of mutations of rearrangement signature 3, 4) number of
mutations of rearrangement signature 5, 5) HRD LOH index, 6) number of
mutations of substitution signature 8.
- **```plot_HRDLOH_HRDetect_Contributions(...)```**: uses the
```HRDetect_pipeline``` output to generate a figure with three plots:
the HDR-LOH index for each sample, the HRDetect BRCAness probability
score for each sample, and the contribution of each of the six features
to the HRDetect BRCAness probability score.
- **```plot_HRDetect_overall```**: uses the ```HRDetect_pipeline```
output to generate an overall plot of the HRDetect BRCAness probability
score.
- **```plot_HRDetect_BootstrapScores```**: overall plot of scores obtained
from the ```HRDetect_pipeline``` output when HRDetect with bootstrap is
enabled.

Function for data visualisation:

- **```genomePlot(...)```**: generates a plot for the visualisation of
somatic variants across the genome, organised in a circle. Variants
plotted are single nucleotide variations (SNV), small insertions and
deletions (indels), copy number variations (CNV) and rearrangements.
- **```plotSignatures(...)```**: this function will plot signatures or
catalogues trying to identify the appropriate mutation type (SNV, DNV, SV...)
from the input row names.

<a name="scripts"/>

## Command line scripts

We provide command line scripts, which can
be used instead of writing your own R code. You can find these scripts in
the ```scripts``` directory.

Currently available scripts are:

- **signatureFit**: mutational signatures analysis using Fit or FitMS. This is a wrapper for the ```signatureFit_pipeline``` R function.
- **hrDetect**: HRDetect pipeline script. This is a wrapper for the ```HRDetect_pipeline``` R function.

You can find user manuals for these command line scripts with detailed
explanation of parameters and examples in the ```userManuals``` folder.

<a name="examples"/>

## Examples

<a name="examplestests"/>

### Test examples

A good place to look for examples is the tests/testthat/ directory,
where for each function of the package we provide a test using test
data. These are the tests that are run when running:

```
devtools::test()
```

Moreover, examples of typical workflows are given below.

<a name="examplese01"/>

### Example 01

In this example we illustrate a typical workflow for signature fit and
HRDetect using two test samples. This code can be easily adapted to be
used on your data, provided you have formatted the input data as
described in the ```?function``` documentation.

This example, along with expected output files, can be found in the
```examples/Example01/``` directory.

We begin by setting variables containing file names. These should be
vectors indexed by sample names.

```
#set directory to this file location
#
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
```

We can now load the SNV tab data and build the SNV mutational
catalogues.

```
#load SNV data and convert to SNV mutational catalogues
SNVcat_list <- list()
for (i in 1:length(SNV_tab_files)){
  tmpSNVtab <- read.table(SNV_tab_files[i],sep = "\t",
                          header = TRUE,check.names = FALSE,
                          stringsAsFactors = FALSE)
  #convert to SNV catalogue, see ?tabToSNVcatalogue or
  #?vcfToSNVcatalogue for details
  res <- tabToSNVcatalogue(subs = tmpSNVtab,genome.v = "hg19")
  colnames(res$catalogue) <- sample_names[i]
  SNVcat_list[[i]] <- res$catalogue
}
#bind the catalogues in one table
SNV_catalogues <- do.call(cbind,SNVcat_list)
```

At this point you can plot the mutational catalogues and compare them
with the expected output files in the ```Example01``` directory.

```
#the catalogues can be plotted as follows
plotSubsSignatures(signature_data_matrix = SNV_catalogues,
                   plot_sum = TRUE,output_file = "SNV_catalogues.pdf")

```

We can now perform signature fit analysis, i.e. identify which
mutational signatures are present in the sample. We suggest to use a
small *a priori* set of signatures, for example here we use the
signatures identified in breast cancer.

```
#fit the 12 breast cancer signatures using the bootstrap signature fit approach
sigsToUse <- c(1,2,3,5,6,8,13,17,18,20,26,30)
subs_fit_res <- Fit(catalogues = SNV_catalogues,
                    signatures = COSMIC30_subs_signatures[,sigsToUse],
                    useBootstrap = TRUE,
                    nboot = 100,
                    nparallel = 4)
plotFit(subs_fit_res,outdir = "signatureFit/")

#The signature exposures can be found here and correspond to the median
#of the bootstrapped runs followed by false positive filters. See
#?Fit for details
snv_exp <- subs_fit_res$exposures
```
In this case we used the ```Fit```
function with ```useBootstrap=TRUE```, which uses the bootstrap fitting method,
and the ```plotFit``` function, which provides
several plots that can be compared to the expected output plots in the
```Example01``` directory.

Finally, we can apply HRDetect on these two samples. Notice that The
HRDetect pipeline allows us to specify the amount of SNV Signature 3 and
8 that we have already estimated above, so we only need to supply the
additional files to compute indels classification, rearrangements and
the copy number based score HRD-LOH. Please see the documentation in
```?HRDetect_pipeline``` for more details.

```
#The HRDetect pipeline will compute the HRDetect probability score for
#the samples to be Homologous Recombination Deficient. HRDetect is a
#logistic regression classifier that requires 6 features to compute the
#probability score. These features can be supplied directly in an input
#matrix, or pipeline can compute these features for you if you supply
#the file names. It is possible to supply a mix of features and file
#names.

#Initialise feature matrix
col_hrdetect <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
input_matrix <- matrix(NA,nrow = length(sample_names),
                       ncol = length(col_hrdetect),
                       dimnames = list(sample_names,col_hrdetect))

#We have already quantified the amount of SNV signatures in the samples,
#so we can supply these via the input matrix
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

```

Also in this case, you can compare your output with the expected output
in the ```Example01``` directory.

<a name="examplese02"/>

### Example 02

In this example we illustrate a typical workflow for signature fit
using organ-specific signatures. This code can be easily adapted to be
used on your data, provided you have formatted the input data as
described in the ```?function``` documentation.

This example, along with expected output files, can be found in the
```examples/Example02/``` directory.

We begin by setting variables containing file names. These should be
vectors indexed by sample names.

```
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
```

We can now load the SNV tab data and build the SNV mutational
catalogues.

```
#load SNV data and convert to SNV mutational catalogues
SNVcat_list <- list()
for (i in 1:length(SNV_tab_files)){
  tmpSNVtab <- read.table(SNV_tab_files[i],sep = "\t",header = TRUE,
                          check.names = FALSE,stringsAsFactors = FALSE)
  #convert to SNV catalogue, see ?tabToSNVcatalogue or ?vcfToSNVcatalogue for details
  res <- tabToSNVcatalogue(subs = tmpSNVtab,genome.v = "hg19")
  colnames(res$catalogue) <- sample_names[i]
  SNVcat_list[[i]] <- res$catalogue
}
#bind the catalogues in one table
SNV_catalogues <- do.call(cbind,SNVcat_list)
```

Then we can plot the catalogues, in this case using pdf format.

```
#the catalogues can be plotted as follows
plotSubsSignatures(signature_data_matrix = SNV_catalogues,
                   plot_sum = TRUE,output_file = "SNV_catalogues.pdf")
```

Than we can perform signature fit using organ specific signatures, here breast cancer signatures.
Typing ```?getOrganSignatures``` will show more information, for example which organs are available.

```
#fit the organ-specific breast cancer signatures using the bootstrap signature fit approach
sigsToUse <- getOrganSignatures("Breast",typemut = "subs")
subs_fit_res <- Fit(catalogues = SNV_catalogues,
                    signatures = sigsToUse,
                    useBootstrap = TRUE,
                    nboot = 100,
                    nparallel = 4)
plotFit(subs_fit_res,outdir = "signatureFit/")

#The signature exposures can be found here and correspond to the median of the bootstrapped
#runs followed by false positive filters. See ?Fit for details
snv_exp <- subs_fit_res$exposures
```

We can now convert the organ-specific signature exposures into reference signature exposures.
This allows us to compare results across different organs, and to perform additional analysis
such as HRDetect, where RefSig 3 and RefSig 8 exposures can be used in place of COSMIC signatures 3 and 8 exposures.

```
#Convert the organ-specific signature exposures into reference signature exposures
snv_exp <- convertExposuresFromOrganToRefSigs(expMatrix = t(snv_exp[,1:(ncol(snv_exp)-1)]),typemut = "subs")

#write the results
writeTable(snv_exp,"RefSigSubsExposures.tsv")
```

<a name="examplese03"/>

### Example 03

In this example we show how to use the multi-step signature fit function along with the Gini-based exposure filter.

Similarly to Examples 01 and 02, we construct the mutational catalogues and then we just need to specify the organ of interest:

```
#perform signature fit using a multi-step approach where organ-specific
#common and rare signatures are used
subs_fit_res <- FitMS(catalogues = SNV_catalogues,
                      exposureFilterType = "giniScaledThreshold",
                      useBootstrap = TRUE,
                      organ = "Breast")
plotFitMS(subs_fit_res,outdir = "signatureFit/")
```

The function FitMS will select automatically the common and rare signatures to use according to the specified organ. As a first step it will fit the common signatures,
and as a second step it will attempt to determine the presence of the rare signatures.

In this example we have also specified to use bootstrap and to use the giniScaledThreshold exposure filter method.
Results can be plotted with the ```plotFitMS``` function.

A manual for FitMS can be found in the ```userManuals``` folder.

<a name="examplese04"/>

### Example 04

Here we provide an example of using the ```signatureFit_pipeline``` function introduced in v2.1.0. This function is meant to automate
some recurrent tasks in signature analysis, such as catalogues generation and selection of signatures to fit.

Similarly to the examples above, we store the location of the files containing the single nucleotide variants into the ```SNV_tab_files``` variable:

```
#set sample names
sample_names <- c("sample1","sample2")
#set the file names.
SNV_tab_files <- c("../../tests/testthat/test_hrdetect_1/test_hrdetect_1.snv.simple.txt",
                   "../../tests/testthat/test_hrdetect_2/test_hrdetect_2.snv.simple.txt")
#name the vectors entries with the sample names
names(SNV_tab_files) <- sample_names
```

Then, we simply call the ```signatureFit_pipeline``` function:

```
pipeline_subs_res <- signatureFit_pipeline(SNV_tab_files = SNV_tab_files,
                                           organ = "Breast",genome.v = "hg19",
                                           fit_method = "FitMS",nparallel = 2)
```

In one function call, the SNV catalogues have been generated, and FitMS has been applied. We can now save the
catalogues, as well as saving the annotated mutations and plotting the results from FitMS:

```
#save catalogues plots
plotSignatures(pipeline_subs_res$catalogues,
               output_file = "SNV_catalogues.pdf",
               ncolumns = 2)
#write annotated SNVs to file
writeTable(pipeline_subs_res$annotated_mutations,"annotated_SNVs.tsv")
#plot the FitMS results
plotFitResults(pipeline_subs_res$fitResults,outdir = "signatureFit/")
```

The function ```plotSignatures``` will infer the type of mutations and select the correct signatures plot function.
In this case, ```plotSignatures``` will notice that these are SNV catalogues and use the ```plotSubsSignatures```
function. Moreover, the function ```plotFitResults``` will infer the fit method (Fit or FitMS) and use the
appropriate plot function, in this case ```plotFitMS```.

The ```signatureFit``` script that can be found in the ```scripts``` folder is a command line wrapper for the
```signatureFit_pipeline``` function.


## Frequently Asked Questions

<a name="faq"/>

**Can I use your package on my Whole Exome Sequenced data?**

Mostly, no. All algorithms included in this package, from signature extraction, to signature fit, to HRDetect, have been developed and tested using only Whole Genome Sequenced data, and are meant to be used only on Whole Genome Sequenced data.
These are some of the reasons why:

- WES have 100x less mutations per sample than WGS
- Using our bootstrap based fitting on WES is not ideal, because bootstrapping WES data can change the profile a lot
- Fitting WGS signatures into WES samples, where the frequency of the trinucleotide contexts are different can also affect the results
- You can convert the WGS signatures into WES by adjusting for trinucleotide context frequency, but the result is far from ideal. It would be better to have signatures extracted from WES to fit into WES samples
- If you then plan to extract WES signatures, this may also be challenging possibly avoiding bootstrapping here as well, and most likely requiring a very large number of WES samples, compared to how many WGS samples are needed to extract WGS signatures
- In the case of HRDetect, the parameters of the classifier have been tuned to Whole Genome data, so it would require retraining using WES data. We did that, and while the performance of a WES HRD classifier was not too bad on breast cancer, we simply discourage to use this for other tissue types. It is just bound to be very unreliable, because the classification would be supported mainly by signature 3 exposures and the number of deletions with microhomology. Fitting signature 3 to WES can be very unreliable, especially in non breast cancer samples, and also the number of indels is usually very low for each WES sample, so one false positive or false negative deletion can make the difference. On top of this, without a standardised data quality control, mutation calling and filtering, there are just too many variables to make this worthwhile, in our opinion. On the other hand, we have plenty of evidence that HRDetect, trained on WGS data using six different mutational signatures, is robust and reliable. 

This said, you may still achieve something decent if you limit yourself to breast cancer, or other organs where we know fairly well what signatures to expect. And you should be able to spot most MMR samples and other hypermutated samples, simply because they have a much higher number of SNVs and Indels.


**In R/HRDetect.R, the function applyHRDetectDavies2017 seems to have mean/sd values hardcoded for each input feature. This results in input data not being standardised to mean=0/sd=1. Is this intentional?
I'm assuming these values come from the results for the 371 samples in the 2017 paper?**

Yes it is intentional. In order for the parameters of the linear model to be meaningful the new input data has to be processed in the same way as the training data, which means that the same mean and sd of the training data should be provided. Please check the R documentation of the function for the details of how to input your data.

**We have a set of samples from non-BRCA cancers, sequenced on a Novaseq, and aligned/called through our in-house pipeline. Are these differences in sequencing/processing large enough to render the parameters of the HRDetect model meaningless?**

Before running HRDetect we make sure that the data is of a good quality and try to reduce artefacts and other false positives. We mostly worked with Sanger cgp pipeline, and also adapted output from other pipelines sometimes with filters to get the same level of specificity, and found that HRDetect is quite robust.
What I mean is that rather than retraining HRDetect for different pipelines, we usually prefer to work on the pipeline until we are happy with the calls.
I can advise to have a look at samples in the same tumour type from highly curated resources like PCAWG. You can browse some at https://signal.mutationalsignatures.com.

**Should I use COSMIC signatures or tissue-specific signatures with my data? Which will give more accurate signature assignment?**

In general, we expect tissue-specific signatures to be more accurate, as we have also demonstrated using simulations in Degasperi et al. 2022, *Science*, Figure S53, though in practice it may depend on the tumour type. Some of the tumour types have more reliable signatures than others (e.g. due to sample size, number of signatures present, similarity between signatures present...). Please double check on SIGNAL that you are happy with the signatures of your cancer type. From version 2.0 of signature.tools.lib, the signature fit algorithm ```FitMS``` automatically selects suitable common and rare signatures for fitting if the parameter ```organ``` is used. From version 2.1.0, the ```signatureFit_pipeline``` also allows the automatic selection of signatures to fit based on tumour types, and also allows to request COSMIC signatures. If COSMIC signatures are used, we advise to select a subset of the COSMIC signatures that are believed to be present in the tumour type to be analysed, rather than attempting to fit all signatures at once.

**How do I interpret the Overall Metrics plot that I obtain from SignatureExtraction? What are all these metrics?**

During the investigation of signature extraction methods, we have defined multiple useful metrics. Some of these were not reported in Degasperi et al. 2020, *Nature Cancer*, as we did not find them essential. The metrics are:

- **norm.Error**: this is the normalised average error between the model and the bootstrapped catalogue. This error is always decreasing as the number of signatures extracted increases.
- **norm.Error (orig. cat.)**: this is the normalised average error between the model and the original catalogue. This error will eventually increase when the number of signatures used is too large, as the algorithm overfits the bootstrapped catalogues.
- **Ave.SilWid**: this is the average silhouette width, which indicates a good the clustering is. The value will be affected by the clustering algorithm used. The Ave.SilWid will begin decreasing when highly similar signatures are extracted, possibly indicating that one is trying to use more signatures than there are actually present in the data.
- **mmcs**: minimum medoid cosine similarity. This indicates the cosine similarity between the two most similar medoids of the clusters. This is in fact the similarity of the two most similar signatures. 
- **min.Min.WCCS**: minimum of the minimum within cluster cosine similarity. This is a goodness of clustering metric. We compute the cosine similarity between the members of each cluster, and for each cluster we report the minimum cosine similarity. Finally we report the minimum of the minimum cosine similarities, which indicates how different the signatures can be in the least homogeneous (most spread) cluster. This value becomes low when at least one cluster is beginning to include very different signatures, possibly indicating complex solutions (see PCA rings in Extended Data Figure 1, in Degasperi et al. 2020, *Nature Cancer*).
- **max.Max.BCCS**: maximum of the maximum between cluster cosine similarity. This is a goodness of clustering metric. We compute the cosine similarity between the members of different clusters, and for each cluster pairs we report the maximum cosine similarity. Finally we report the maximum of the maximum cosine similarities, which indicates how similar the signatures can be when they belong to different clusters. This value becomes high when two clusters are beginning to include very similar signatures, possibly indicating complex solutions (see PCA rings in Extended Data Figure 1, in Degasperi et al. 2020, *Nature Cancer*).

While it can be interesting to look at all these metrics, we recommend to mostly use **norm.Error (orig. cat.)** and **Ave.SilWid**, using a suitable number of bootstraps and repeats and with the Clustering with Matching (MC) clustering algorithm, as described in Degasperi et al. 2020, *Nature Cancer*.

**What parameters should I use for signature extraction?**

It depends on how many samples you have, with more samples requiring more repeats. In general, we advise to use 20 bootstraps, 200 repeats, the clustering with matching algorithm (CM), the KLD objective function (nmfmethod brunet) and RTOL=0.001. The number of repeats (nrepeats) is more important than the number of bootstraps, because at least 100 repeats are necessary to obtain enough results to then select best runs using the RTOL threshold (Extended Data Figure 1 of Degasperi et al. 2020, *Nature Cancer*). The disadvantage of this method is that this can increase the computation time required significantly if one uses more than 500/1000 samples. 

**Where do I find the activity (exposures) matrix after running the signature extraction?**

We do not report the activity matrix from the extraction, only the
signatures. We then use a signature fit procedure (here the function ```Fit``` or ```FitMS```)
to estimate the activities. We tend to think of the estimation of the
activities as a separate procedure, performed holding a given set of a
priori signatures fixed. This allows for more playing around with the
signature fit, for example excluding clear false positive signatures from
the fit.

**Which Structural Variants caller and filters should I use?**

Our team has been using mostly Manta and Brass. For Manta, we tend to tune false positives using the PR column.

https://github.com/Illumina/manta/tree/master/docs/userGuide

After running Manta you will have to convert your output to bedpe using a vcf to bedpe converter. There are a few out there, you can try this one:

https://github.com/arq5x/lumpy-sv/blob/master/scripts/vcfToBedpe



