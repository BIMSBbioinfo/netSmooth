<div align="center">
	<img src="hex-netsmooth.png" alt="netsmooth"/>
</div>

netSmooth
---------

[![Build Status](https://travis-ci.org/BIMSBbioinfo/netSmooth.svg?branch=master)](https://travis-ci.org/BIMSBbioinfo/netSmooth) [![codecov](https://codecov.io/gh/BIMSBbioinfo/netSmooth/branch/master/graph/badge.svg)](https://codecov.io/gh/BIMSBbioinfo/netSmooth)

_netSmooth: A Network smoothing based method for Single Cell RNAseq and other gene-based genomics data sets_
netSmooth is an R package for network smoothing of single cell RNA sequencing data. Using gene interaction networks such as protein-
protein interactions as priors for gene co-expression, netsmooth improves cell type identification from noisy, sparse scRNAseq data.
The smoothing method is suitable for other gene-based omics data sets as well such as proteomics, copy-number variation, etc.

### Installation

	library(devtools)
	install_github("BIMSBbioinfo/netSmooth")

### Usage
For detailed usage information see  [the vignette](http://htmlpreview.github.io/?https://github.com/BIMSBbioinfo/netSmooth/blob/master/vignettes/netSmoothIntro.html). In addition,
R package has full function documentation with examples. 

### Contributing

Fork and send a pull request. Or just e-mail us.

-------------------------
@jonathanronen, BIMSBbioinfo, 2017

