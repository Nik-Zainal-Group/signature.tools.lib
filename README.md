# Signature Tools Lib R package

## Table of content

- [Introduction to the package](#intro)
- [Systems Requirements](#req)
- [Testing the package](#test)
- [How to use this package](#howtouse)
- [Package documentation](#docs)
- [Functions provided by the package](#functions)
- [Examples](#examples)
		- [Test examples](#examplestests)
		- [Example 01](#examplese01)


<a name="intro"/>

## Introduction to the package

Signature Tools Lib is an R package for mutational signatures analysis. Mutational signatures are patterns of somatic mutations that reveal what mutational processes have been active in a cell. Mutational processes can be due to exposure to mutagens, such as chemicals present in cigarettes, or defects in DNA repair pathways, such as Homologous Recombination Repair.

The package supports hg19 and hg38 as well as mm10. It provides our latest algorithms for signature fit and extraction, as well as various utility functions and the HRDetect pipeline. The list and description of these functions is given below.

<a name="req"/>

## Systems Requirements

No special hardware is required to run this software. This is an R package so it will work on any computer with R (>=3.2.1) installed. You can install this package by entering the R environment from the main directory and typing:

```
install.packages("devtools")
devtools::install()
```

The above commands will take care of all the dependencies. In some cases, it may happen that some dependencies cannot be installed automatically. In this case, an error will indicate which package could not be found, then simply install these packages manually. The installation should not take more than a few minutes.

This is the full list of R package dependencies:

```
    VariantAnnotation,
    BSgenome.Hsapiens.UCSC.hg38,
    BSgenome.Hsapiens.1000genomes.hs37d5,
    BSgenome.Mmusculus.UCSC.mm10,
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

You can test the package by entering the package main directory and typing from the R invironment:

```
devtools::test()
```

<a name="howtouse"/>

## How to use this package

**PLEASE NOTE:** project-specific file conversions and filtering should be done before and outside the use of this library. This library should only report to the user the data format errors and suggest how to correct them, while it is the responsibility of the user to supply the necessary data formatted correctly with the necessary features and correct column names.

<a name="docs"/>

## Package documentation

**DOCUMENTATION:** Documentation for each of the functions below is provided as R documentation, and it is installed along with the R package. The documentation should give detailed explanation of the input data required, such as a list of data frame columns and their explanation. To access the documentation you can use the ```?function``` syntax in R, for each of the functions below. For example, in R or RStudio, type ```?HRDetect_pipeline```.

<a name="functions"/>

## Functions provided by the package

Functions for file conversion/manipulation:

- **```bedpeToRearrCatalogue(...)```**: converts a data frame (loaded from a BEDPE file) into a rearrangement catalogue. If the columns ```is.clustered``` or ```svclass``` are not present, this function will attempt to compute them.
- **```vcfToSNVcatalogue(...)```**: converts a VCF file into a SNV catalogue.
- **```tabToSNVcatalogue(...)```**: converts a data frame into a SNV catalogue. The data frame should have the following minimal columns: chr, position, REF, ALT.

Functions for signature extraction and signature fit

- **```SignatureFit_withBootstrap(...)```**: fit a given set of mutational signatures into mutational catalogues to extimate the activty/exposure of each of the given signatures in the catalogues. Implementation of method similar to Huang 2017, Detecting presence of mutational signatures with confidence, which uses a bootstrap apporach to calculate the empirical probability of an exposure to be larger or equal to a given threshold (i.e. 5% of mutations of a sample). This probability can be used to decide which exposures to remove from the initial fit, thus increasing the sparsity of the exposures.
- **```SignatureFit_withBootstrap_Analysis(...)```**: this function is a wrapper for the function ```SignatureFit_withBootstrap```, which produces several plots for each sample in the catalogues
- **```SignatureExtraction(...)```**: perform signature extraction, by applying NMF to the input matrix. Multiple NMF runs and bootstrapping is used for robustness, followed by clustering of the solutions. A range of number of signatures to be used is required.
- **```plotSubsSignatures(...)```**: function to plot one or more substitution signatures or catalogues.
- **```plotRearrSignatures(...)```**: function to plot one or more rearrangement signatures or catalogues.

Functions for HRD indexes

- **```ascatToHRDLOH(...)```**: compute the HRD-LOH index from a data frame formatted similarly to an ASCAT output file. This is a wrapper for the function ```calc.hrd```, with a simplified interface, where only a data frame is requested.
- **```calc.hrd(...)```**: compute HRD-LOH (Loss of Heterozygosity), written by Nicolai Juul Birkbak
- **```calc.ai(...)```**: compute HRD-NtAI (Number of telomeric Allelic Imbalances), written by Nicolai Juul Birkbak
- **```calc.lst(...)```**: compute HRD-LST (Large-scale state transitions), written by Nicolai Juul Birkbak

Functions for Indels Classification

- **```vcfToIndelsClassification(...)```**: converts a VCF file containing indels into a data frame where the indels are classified. Also returns a summary of count and proportion of the various classes of indels.
- **```tabToIndelsClassification(...)```**: converts a data frame containing indels into a data frame where the indels are classified. Also returns a summary of count and proportion of the various classes of indels.

Functions for HRDetect

- **```HRDetect_pipeline(...)```**: this is a flexible interface to the HRDetect pipeline to compute the HRDetect BRCAness probability score, as published in Davies et al. 2017.
- **```applyHRDetectDavies2017(...)```**: given a data frame with samples as rows and features as columns, this function will compute the HRDetect BRCAness probability for each sample. The following six features should be included in the matrix: 1) proportion of deletions with micro-homology, 2) number of mutations of substitution signature 3, 3) number of mutations of rearrangemet signature 3, 4) number of mutations of rearrangemet signature 5, 5) HRD LOH index, 6) number of mutations of substitution signature 8.
- **```plot_HRDLOH_HRDetect_Contributions(...)```**: uses the ```HRDetect_pipeline``` output to generate a figure with three plots: the HDR-LOH index for each sample, the HRDetect BRCAness probability score for each sample, and the contribution of each of the six features to the HRDetect BRCAness probability score.
- **```plot_HRDetect_overall```**: uses the ```HRDetect_pipeline``` output to generate an overall plot of the HRDetect BRCAness probability score.

Function for data visualisation

- **```genomePlot(...)```**: generates a plot for the visualisation of somatic variants across the genome, organised in a circle. Variants plotted are single nucleotide variations (SNV), small insertions and deletions (indels), copy number variations (CNV) and rearrangements.

Function for web formats export

- **```export_SignatureFit_withBootstrap_to_JSON```**: Given a res file obtained from the ```SignatureFit_withBootstrap``` or ```SignatureFit_withBootstrap_Analysis``` function, export it to a set of JSON files that can be used for web visualisation

<a name="examples"/>

## Examples

<a name="examplestests"/>

### Test examples

A good place to look for examples is the tests/testthat/ directory, where for each function of the package we provide a test using test data. These are the tests that are run when running:

```
devtools::test()
```

Moreover, examples of typical workflows are given below.

<a name="examplese01"/>

### Example 01

In this example we illustrate a typical workflow for signature fit and HRDetect using two test samples. This code can be easily adapted to be used on your data, provided you have formatted the input data as described in the ```?function``` documentation.

This example, along with expected output files, can be found in the ```examples/Example01/``` directory.

We begin by setting variables containing file names. These should be vectors indexed by sample names.

```
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
```

