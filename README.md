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

