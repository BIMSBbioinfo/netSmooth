<div align="center">
	<img src="hex-netsmooth.png" alt="netsmooth"/>
</div>


---------

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1119064.svg)](https://doi.org/10.5281/zenodo.1119064)
[![Build Status](http://bioconductor.org/shields/build/devel/bioc/netSmooth.svg)](https://travis-ci.org/BIMSBbioinfo/netSmooth) [![codecov](https://codecov.io/gh/BIMSBbioinfo/netSmooth/branch/master/graph/badge.svg)](https://codecov.io/gh/BIMSBbioinfo/netSmooth) [![BioC_years](http://www.bioconductor.org/shields/years-in-bioc/netSmooth.svg)](http://www.bioconductor.org/packages/release/bioc/html/netSmooth.html)

**_netSmooth_: A Network smoothing based method for single cell RNA-seq**
-----
_netSmooth_ is an R package for network smoothing of single cell RNA sequencing data. Using gene interaction networks such as protein-
protein interactions as priors for gene co-expression, _netsmooth_ improves cell type identification from noisy, sparse scRNA-seq data.
The smoothing method is suitable for other gene-based omics data sets such as proteomics, copy-number variation, etc.

The algorithm uses a network-diffusion based approach which takes in a network (such as PPI network) and gene-expression matrix. The gene 
expression values in the matrix are smoothed using the interaction information in the network. The network-smoothing parameter is 
optimized using a robust clustering approach.

For a detailed exposition, check out [our paper on F1000Research](https://f1000research.com/articles/7-8/v2).

### Installation

_netSmooth_ is available via Bioconductor:

	if (!requireNamespace("BiocManager", quietly=TRUE))
    	install.packages("BiocManager")
	BiocManager::install("netSmooth")

Alternatively, using `devtools`:

	library(devtools)
	install_github("BIMSBbioinfo/netSmooth")

### Usage
For detailed usage information see  [the vignette](http://htmlpreview.github.io/?https://github.com/BIMSBbioinfo/netSmooth/blob/master/vignettes/netSmoothIntro.html). In addition,
the R package has full function documentation with examples.

### How to cite
Please cite the _netSmooth_ paper:

> Ronen J and Akalin A. _netSmooth_: Network-smoothing based imputation for single cell RNA-seq [version 2; referees: 2 approved]. F1000Research 2018, 7:8 (doi: 10.12688/f1000research.13511.2)

### License

_netSmooth_ is available under a GPLv3 license.

### Contributing

Fork and send a pull request. Or just e-mail us.

-------------------------
@jonathanronen, BIMSBbioinfo, 2017

