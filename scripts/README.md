# signature.tools.lib scripts

## Introduction

This folder includes scripts to perform analysis using the R package ```signature.tools.lib```
in unix command line style, without having to write R code.

Currently available scripts are:

- **signatureFit**: mutational signatures analysis using Fit or FitMS. This is a wrapper for the ```signatureFit_pipeline``` R function.
- **hrDetect**: HRDetect pipeline script. This is a wrapper for the ```HRDetect_pipeline``` R function.
- **genomeChart**: visualisation of somatic variants using a circle diagram and other graphs like mutational catalogues.
- **solutionSelectionForFitMS**: change selection criteria for FitMS alternative rare signature solutions and/or select solutions manually.

You can use the flag ```-h``` to access a list of options for each script, for example ```hrDetect -h```.

### Install in your PATH

You can move these scripts to a directory in your PATH, so that it can be used from
the command line like any other program.

For example, you can create a ```bin``` directory in your home directory:

```
mkdir ~/bin
```

And add it to your path:

```
export PATH=~/bin:$PATH
```

You can add the above line to your ```~/.bashrc``` file so you don't have to type it every time.
Then, you can simply copy the files:

```
cd /path/to/where/you/downloaded/signature.tools.lib/scripts/
cp signatureFit ~/bin/
cp hrDetect ~/bin/
cp genomeChart ~/bin/
cp solutionSelectionForFitMS ~/bin/
```

Alternatively, you can create a symbolic link:

```
cd ~/bin/
ln -s /path/to/where/you/downloaded/signature.tools.lib/scripts/signatureFit
ln -s /path/to/where/you/downloaded/signature.tools.lib/scripts/hrDetect
ln -s /path/to/where/you/downloaded/signature.tools.lib/scripts/genomeChart
ln -s /path/to/where/you/downloaded/signature.tools.lib/scripts/solutionSelectionForFitMS
```

You should now have the scripts working anywhere in your command line.

