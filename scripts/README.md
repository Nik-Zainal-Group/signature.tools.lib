# Signature Tools Lib Scripts

## Introduction

This folder includes scripts to perform analysis using the R package ```signature.tools.lib```
in unix command line style, without having to write R code.

## signatureFit

For the help on how to use this script type:

```
./signatureFit.R -h
```

## hrDetect

For the help on how to use this script type:

```
./hrDetect.R -h
```

### Install in your PATH

You can move this script to a directory in your PATH, so that it can be used from
the command line like any other program.

For example, you can create a ```bin``` directory in your home directory:

```
mkdir ~/bin
```

And add it to your path:

```
export PATH=~/bin:$PATH
```

You can add the above line to your ```~/.bashrc``` file so you don't have to retype it every time.
Then, you can simply copy the files and even remove the ```.R``` extension:

```
cp signatureFit.R ~/bin/signatureFit
cp hrDetect.R ~/bin/hrDetect
```

You should now have the scripts working anywhere in your command line.
