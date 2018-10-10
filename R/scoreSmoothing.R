#' Calculate a score for a smoothing result, for picking the best alpha value
#' @param x    the network-smoothed expression matrix
#' @param method    the scoring method. 'entropy' calculates shannon entropy
#'                  in a 2D PCA of the data. 'robustness' performs robsut
#'                  clustering and reports the proportion of samples in
#'                  robust clusters
#' @return the score
#' @keywords internal
scoreSmoothing <- function(x, method=c('entropy', 'robustness'),
    is.counts=TRUE, ...) {
    if(class(x)=='matrix' || any(is(x)=='Matrix' || class(x)=='DelayedMatrix')) {
        x <- x
        se <- SummarizedExperiment::SummarizedExperiment(x)
    } else if(class(x)=='SummarizedExperiment') {
        se <- x
        x <- assay(x)
    } else stop("must be matrix/Matrix or SummarizedExperiment object")

    method <- match.arg(method)
    if(method=='entropy') {
        score <- calc2DEntropy(dimReduce(x,
            flavor='pca', is.counts=is.counts))
    }
    else if(method=='robustness') {
        score <- robustClusters(se, runMergeClusters=FALSE,
            ...)$proportion.robust
    }
    return(score)
}
