#' Get lower dimension embedding
#'
#' @param expr    gene expresison matrix [GENES x SAMPLES]
#' @param flavor    the algorithm to use to obtain the dimensionality reduction
#'                  must be in c('pca', 'tsne')
#' @param k    the number of dimensions in the reduced dimension representation
#' @param is.counts    logical: is `expr` counts data
#' @return    reduced dimensionality representation
dimReduce <- function(expr, flavor=c('pca', 'tsne'), k=2, is.counts=TRUE) {
    flavor <- match.arg(flavor)
    if(flavor=='pca') function.to.call <- scater::plotPCA
    if(flavor=='tsne') function.to.call <- scater::plotTSNE

    if(is.counts) sce <- scater::newSCESet(countData = expr)
    else sce <- suppressWarnings( scater::newSCESet(exprsData = expr) )

    sce <- function.to.call(sce, return_SCESet=TRUE, draw_plot=FALSE,
                            ncomponents=k)
    red.dim <- data.frame(reducedDimension(sce))[,1:k]
    colnames(red.dim) <- paste0(flavor, 1:k)
    return(red.dim)
}

#' Calculate entropy in 2D data
#' @param x    the 2D data to get entropy from
#' @param numBins1    the number of bins along the first dimension to
#'                    discretize data into
#' @param numBins2    the number of bins along the second dimension to
#'                    discretize data into
#' @return The Shannon entropy in the 2D data x
calc2DEntropy <- function(x, numBins1=20, numBins2=20) {
    embedding.discretized <- entropy::discretize2d(x[,1], x[,2], numBins1, numBins2)
    return(entropy::entropy(embedding.discretized))
}

#' Pick the dimensionality reduction method for a dataset that gives the
#' 2D embedding with the highest entropy
#'
#' @param expr  expression matrix [GENES x SAMPLES]
#' @param flavors    list of dimensionality reduction algorithms to try
#' @param is.counts    logical: is exprs count data
#' @usage pickDimReduction(expr)
#' @return    name of dimensionality reduction method that gives the highest
#'            2d entropy
#' @export
pickDimReduction <- function(expr, flavors=c('pca', 'tsne'), is.counts=TRUE) {
    entropies <- sapply(flavors, function(flavor) {
        calc2DEntropy(dimReduce(expr, flavor=flavor, is.counts=is.counts))
    })
    return(names(which.max(entropies)))
}
