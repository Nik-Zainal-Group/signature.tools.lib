# Signature Tools Lib R package

You can install this R package by entering the main directory and typing:

```
install.packages("devtools")
devtools::install()
```

You can also test the package by typing:

```
devtools::test()
```

Useful commands for development workflow can be found in the file ```signatures.tools.lib.develop.R```.


List of functions:

- bedpeToRearrCatalogue(): converts a data frame (loaded from a BEDPE file) into a rearrangement catalogue. If the columns ```is.clustered``` or ```svclass``` are not present, this function will attempt to compute them.
- vcfToSNVcatalogue(): converts a VCF file into a SNV catalogue.
- SignatureFit_withBootstrap(): Fit a given set of mutational signatures into mutational catalogues to extimate the activty/exposure of each of the given signatures in the catalogues. Implementation of method similar to Huang 2017, Detecting presence of mutational signatures with confidence, which uses a bootstrap apporach to calculate the empirical probability of an exposure to be larger or equal to a given threshold (i.e. 5% of mutations of a sample). This probability can be used to decide which exposures to remove from the initial fit, thus increasing the sparsity of the exposures.
- SignatureExtraction(): Perform signature extraction, by applying NMF to the input matrix. Multiple NMF runs and bootstrapping is used for robustness, followed by clustering of the solutions. A range of number of signatures to be used is required.
