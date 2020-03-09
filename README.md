# Signature Tools Lib R package

## Table of content

- [Introduction to the package](#intro)
- [How to cite us](#cite)
- [Systems Requirements](#req)
- [Testing the package](#test)
- [How to use this package](#howtouse)
- [Package documentation](#docs)
- [Functions provided by the package](#functions)
- [Examples](#examples)
  - [Test examples](#examplestests)
  - [Example 01](#examplese01)
- [Frequently Asked Questions](#faq)

<a name="intro"/>

## Introduction to the package

Signature Tools Lib is an R package for mutational signatures analysis.
Mutational signatures are patterns of somatic mutations that reveal what
mutational processes have been active in a cell. Mutational processes
can be due to exposure to mutagens, such as chemicals present in
cigarettes, or defects in DNA repair pathways, such as Homologous
Recombination Repair.

The package supports hg19 and hg38 as well as mm10. It provides our
latest algorithms for signature fit and extraction, as well as various
utility functions and the HRDetect pipeline. The list and description of
these functions is given below.


## How to cite us

<a name="cite"/>

If you are using this package please cite:

A. Degasperi et al. **A practical framework and online tool for mutational signature analyses show intertissue variation and driver dependencies**, *Nature Cancer*, [https://doi.org/10.1038/s43018-020-0027-5], 2020.

In this publication you will find a lot of useful information for using this package.

<a name="req"/>

## Systems Requirements

No special hardware is required to run this software. This is an R
package so it will work on any computer with R (>=3.2.1) installed. You
can install this package by entering the R environment from the main
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
    doMC,
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
    BSgenome
```

<a name="test"/>

## Testing the package

You can test the package by entering the package main directory and
typing from the R invironment:

```
devtools::test()
```

<a name="howtouse"/>

## How to use this package
## 
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

Functions for file conversion/manipulation:

- **```bedpeToRearrCatalogue(...)```**: converts a data frame (loaded
from a BEDPE file) into a rearrangement catalogue. If the columns
```is.clustered``` or ```svclass``` are not present, this function will
attempt to compute them.
- **```vcfToSNVcatalogue(...)```**: converts a VCF file into a SNV
catalogue.
- **```tabToSNVcatalogue(...)```**: converts a data frame into a SNV
catalogue. The data frame should have the following minimal columns:
chr, position, REF, ALT.

Functions for signature extraction and signature fit

- **```SignatureFit_withBootstrap(...)```**: fit a given set of
mutational signatures into mutational catalogues to extimate the
activty/exposure of each of the given signatures in the catalogues.
Implementation of method similar to Huang 2017, Detecting presence of
mutational signatures with confidence, which uses a bootstrap apporach
to calculate the empirical probability of an exposure to be larger or
equal to a given threshold (i.e. 5% of mutations of a sample). This
probability can be used to decide which exposures to remove from the
initial fit, thus increasing the sparsity of the exposures.
- **```SignatureFit_withBootstrap_Analysis(...)```**: this function is a
wrapper for the function ```SignatureFit_withBootstrap```, which
produces several plots for each sample in the catalogues
- **```SignatureExtraction(...)```**: perform signature extraction, by
applying NMF to the input matrix. Multiple NMF runs and bootstrapping is
used for robustness, followed by clustering of the solutions. A range of
number of signatures to be used is required.
- **```plotSubsSignatures(...)```**: function to plot one or more
substitution signatures or catalogues.
- **```plotRearrSignatures(...)```**: function to plot one or more
rearrangement signatures or catalogues.

Functions for HRD indexes

- **```ascatToHRDLOH(...)```**: compute the HRD-LOH index from a data
frame formatted similarly to an ASCAT output file. This is a wrapper for
the function ```calc.hrd```, with a simplified interface, where only a
data frame is requested.
- **```calc.hrd(...)```**: compute HRD-LOH (Loss of Heterozygosity),
written by Nicolai Juul Birkbak
- **```calc.ai(...)```**: compute HRD-NtAI (Number of telomeric Allelic
Imbalances), written by Nicolai Juul Birkbak
- **```calc.lst(...)```**: compute HRD-LST (Large-scale state
transitions), written by Nicolai Juul Birkbak

Functions for Indels Classification

- **```vcfToIndelsClassification(...)```**: converts a VCF file
containing indels into a data frame where the indels are classified.
Also returns a summary of count and proportion of the various classes of
indels.
- **```tabToIndelsClassification(...)```**: converts a data frame
containing indels into a data frame where the indels are classified.
Also returns a summary of count and proportion of the various classes of
indels.

Functions for HRDetect

- **```HRDetect_pipeline(...)```**: this is a flexible interface to the
HRDetect pipeline to compute the HRDetect BRCAness probability score, as
published in Davies et al. 2017.
- **```applyHRDetectDavies2017(...)```**: given a data frame with
samples as rows and features as columns, this function will compute the
HRDetect BRCAness probability for each sample. The following six
features should be included in the matrix: 1) proportion of deletions
with micro-homology, 2) number of mutations of substitution signature 3,
3) number of mutations of rearrangemet signature 3, 4) number of
mutations of rearrangemet signature 5, 5) HRD LOH index, 6) number of
mutations of substitution signature 8.
- **```plot_HRDLOH_HRDetect_Contributions(...)```**: uses the
```HRDetect_pipeline``` output to generate a figure with three plots:
the HDR-LOH index for each sample, the HRDetect BRCAness probability
score for each sample, and the contribution of each of the six features
to the HRDetect BRCAness probability score.
- **```plot_HRDetect_overall```**: uses the ```HRDetect_pipeline```
output to generate an overall plot of the HRDetect BRCAness probability
score.

Function for data visualisation

- **```genomePlot(...)```**: generates a plot for the visualisation of
somatic variants across the genome, organised in a circle. Variants
plotted are single nucleotide variations (SNV), small insertions and
deletions (indels), copy number variations (CNV) and rearrangements.

Function for web formats export

- **```export_SignatureFit_withBootstrap_to_JSON```**: Given a res file
obtained from the ```SignatureFit_withBootstrap``` or
```SignatureFit_withBootstrap_Analysis``` function, export it to a set
of JSON files that can be used for web visualisation

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
                   plot_sum = TRUE,output_file = "SNV_catalogues.jpg")

```

We can now perform signature fit analysis, i.e. identify which
mutational signatures are present in the sample. We suggest to use a
small *a priori* set of signatures, for example here we use the
signatures identified in breast cancer.

```
#fit the 12 breast cancer signatures using the bootstrap signature fit approach
sigsToUse <- c(1,2,3,5,6,8,13,17,18,20,26,30)
subs_fit_res <- SignatureFit_withBootstrap_Analysis(outdir = "signatureFit/",
                                    cat = SNV_catalogues,
                                    signature_data_matrix = COSMIC30_subs_signatures[,sigsToUse],
                                    type_of_mutations = "subs",
                                    nboot = 100,nparallel = 4)

#The signature exposures can be found here and correspond to the median
#of the boostrapped runs followed by false positive filters. See
#?SignatureFit_withBootstrap_Analysis for details
snv_exp <- subs_fit_res$E_median_filtered
```
In this case we used the ```SignatureFit_withBootstrap_Analysis```
function, which will use the bootstrap fitting method and provide
several plots, which can be compared to the expected output plots in the
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
input_matrix[colnames(snv_exp),"SNV3"] <- snv_exp["Signature.3",]
input_matrix[colnames(snv_exp),"SNV8"] <- snv_exp["Signature.8",]

#run the HRDetect pipeline, for more information see ?HRDetect_pipeline
res <- HRDetect_pipeline(input_matrix,
                         genome.v = "hg19",
                         SV_bedpe_files = SV_bedpe_files,
                         Indels_vcf_files = Indels_vcf_files,
                         CNV_tab_files = CNV_tab_files,
                         nparallel = 2)

#save HRDetect scores
writeTable(res$hrdetect_output,file = "HRDetect_res.tsv")

```

Also in this case, you can compare your output with the expected output
in the ```Example01``` directory.



## Frequently Asked Questions

<a name="faq"/>

**In R/HRDetect.R, the function applyHRDetectDavies2017 seems to have mean/sd values hardcoded for each input feature. This results in input data not being standardized to mean=0/sd=1. Is this intentional?
I'm assuming these values come from the results for the 371 samples in the 2017 paper?**

Yes it is intentional. In order for the parameters of the linear model to be meaningful the new input data has to be processed in the same way as the training data, which means that the same mean and sd of the training data should be provided. Please check the R documentation of the function for the details of how to input your data.

**We have a set of samples from non-BRCA cancers, sequenced on a Novaseq, and aligned/called through our in-house pipeline. Are these differences in sequencing/processing large enough to render the parameters of the HRDetect model meaningless?**

Before running HRDetect we make sure that the data is of a good quality and try to reduce artefacts and other false positives. We mostly worked with Sanger cgp pipeline, and also adapted output from other pipelines sometimes with filters to get the same level of specificity, and found that HRDetect is quite robust.
What I mean is that rather than retraining HRDetect for different pipelines, we usually prefer to work on the pipeline until we are happy with the calls.
I can advise to have a look at samples in the same tumour type from highly curated resources like PCAWG. You can browse some at https://signal.mutationalsignatures.com.

**Should I use COSMIC signatures or Tissue Specific signatures with my data? Which will give more accurate signature assignment?**

In general we expect tissue specific signatures to be more accurate, though for now in practice it depends on tumour types.
We performed signature extraction in each organ to get the organ specific signatures, which you can also find here https://signal.mutationalsignatures.com/. Some of the tumour types have more reliable signatures than others (e.g. due to sample size, number of signatures present, similarity between signatures present...). Please double check on SIGNAL that you are happy with the signatures (both substitution and rearrarangement) of your cancer type. We plan to update these signatures as we get more samples, hopefully making them more reliable with the time.
If you choose to use the COSMIC signatures, that can also be fine, though we do not advise fitting all 30 (or 60 with the new PCAWG paper) subs signatures at the same time. It is better to try to find which COSMIC signatures to expect in your tumour type. A disadvantage of using the COSMIC option in the signature_type, is that there is no equivalent of COSMIC for rearrangement signatures, so only the rearrangement signatures from the 560 breast cancer paper are used.
In general, a more advanced way to do this, whichever signatures you choose to fit, would be to check whether your fitted model looks like the data (for example using the ```SignatureFit_withBootstrap_Analysis(...)``` function, which generates a detailed plot for each sample). The difference between model and data can give you a clue of other signatures that should be added in the fitting (if any). Also, if you are unsure whether a signature is really there or not, then you can try to refit without that signature and see if the similarity between model and data is reduced significantly and whether a pattern similar to the removed signature appears in the difference plot.

**How do I interpret the Overall Metrics plot that I obtain from SignatureExtraction? What are all these metrics?**

During the investigation of Signature Extraction methods, we have defined multiple useful metrics. Some of these were not reported in the Nature Cancer paper, as we did not find them essential. The metrics are:

- **norm.Error**: this is the normlised average error between the model and the bootstrapped catalogue. This error is always decreasing as the number of signatures extracted increases.
- **norm.Error (orig. cat.)**: this is the normalised average error between the modela and the original catalogue. This error will eventually increase when the number of signatures used is too large, as the algorithm overfits the bootstrapped catalogues.
- **Ave.SilWid**: this is the average silhouette width, which indicates a good the clustering is. The value will be affected by the clustering algorithm used. The Ave.SilWid will begin decreasing when highly similar signatures are extracted, possibly indicating that one is trying to use more signatures than there are actually present in the data.
- **mmcs**: minimum medoid cosine similarity. This indicates the cosine similarity between the two most similar medoids of the clusters. This is in fact the similarity of the two most similar signatures. 
- **min.Min.WCCS**: minimum of the minimum within cluster cosine similarity. This is a goodness of clustering metric. We compute the cosine similarity between the members of each cluster, and for each cluster we report the minimum cosine similarity. Finally we report the minimum of the minimum cosine similarities, which indicates how different the signatures can be in the least homogeneous (most spread) cluster. This value becomes low when at least one cluster is beginning to include very different signatures, possibly indicating complex solutions (see PCA rings in Extended Data Figure 1, in the Nature Cancer paper).
- **max.Max.BCCS**: maximum of the maximum between cluster cosine similarity. This is a goodness of clustering metric. We compute the cosine similarity between the members of different clusters, and for each cluster pairs we report the maximum cosine similarity. Finally we report the maximum of the maximum cosine similarities, which indicates how similar the signatures can be when they belong to different clusters. This value becomes high when two clusters are beginning to include very similar signatures, possibly indicating complex solutions (see PCA rings in Extended Data Figure 1, in the Nature Cancer paper).

While it can be interesting to look at all these metrics, we reccommend to mostly use **norm.Error (orig. cat.)** and **Ave.SilWid**, using a suitable number of bootstraps and repeats and with the Clustering with Matching (MC) clustering algorithm, as described in our Nature cancer paper.

**What parameters should I use for signature extraction?**

It depends on how many samples you have, with more samples requiring more repeats. In general, we advise to use 20 bootstraps, 200 repeats, the clustering with matching algorithm (CM), the KLD objective function (nmfmethod brunet) and RTOL=0.001. The number of repeats (nrepeats) is more important than the number of bootstraps, because at least 100 repeats are necessary to obtain enough results to then select best runs using the RTOL threshold (Extended Data Figure 1 of our Nature Cancer paper). The disadvantage of this method is that this can increase the computation time required significantly if one uses more than 500/1000 samples. 

**Where do I find the activity (exposures) matrix after running the signature extraction?**

We do not report the activity matrix from the extraction, only the
signatures. We then use a signature fit procedure (described in the Nature Cancer paper, here ```SignatureFit_withBootstrap```)
to estimate the activities. We tend to think of the estimation of the
activities as a separate procedure, performed holding a given set of a
priori signatures fixed. This allows for more playing around with the
signature fit, for example excluding clear false positive signatures from
the fit.
