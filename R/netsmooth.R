#' Calculate a score for a smoothing result, for picking the best alpha value
#' @param x    the network-smoothed expression matrix
#' @param method    the scoring method. 'entropy' calculates shannon entropy
#'                  in a 2D PCA of the data. 'robustness' performs robsut
#'                  clustering and reports the proportion of samples in
#'                  robust clusters
#' @return the score
scoreSmoothing <- function(expr, method=c('entropy', 'robustness'),
                           is.counts=TRUE, ...) {
    if(class(expr)=='matrix') {
        x <- expr
        se <- SummarizedExperiment::SummarizedExperiment(expr)
    } else if(class(expr)=='SummarizedExperiment') {
        x <- assay(expr)
        se <- expr
    } else stop("must be matrix or SummarizedExperiment object")

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

#' Perform network smoothing of gene expression or other omics data
#' @param x    matrix or SummarizedExperiment
#' @param adjMatrix    adjacency matrix of gene network to use
#' @param alpha    numeric in [0,1] or 'audo'. if 'auto', the optimal
#'                 value for alpha will be automatically chosen among the values
#'                 specified in `autoAlphaRange`, using the strategy
#'                 specified in `autoAlphaMethod`
#' @param autoAlphaMethod    if 'robustness', pick alpha that gives the
#'                              highest proportion of samples in robust clusters
#'                              if 'entropy', pick alpha that gives highest
#'                              Shannon entropy in 2D PCA embedding
#' @param autoAlphaRange    if `alpha='optimal'`, search these values
#'                             for the best alpha
#' @param autoAlphaDimReduceFlavor    algorithm for dimensionality reduction
#'                                    that will be used to pick the optimal
#'                                    value for alpha. Either the 2D embedding
#'                                    to calculate the Shannon entropy for (if
#'                                    `autoAlphaMethod='entropy'`), or the
#'                                    dimensionality reduction algorithm to be
#'                                    used in robust clustering (if
#'                                    `autoAlphamethod='robustness'`)
#' @param summarizedExperimentAssay    if `x` is a SummarizedExperiment object,
#'                                     the index of the assay to use
#' @param is.counts    logical: is the assay count data
#' @param numcores    number of parallel processes to use for auto alpha choice
#' @param ...    arguments passed on to `robustClusters` if using the robustness
#'               criterion for optimizing alpha
#' @return network-smoothed gene expression matrix or SummarizedExperiment
#'         object
#' @export
netSmooth <- function(x, adjMatrix, alpha='auto',
                      autoAlphaMethod=c('robustness', 'entropy'),
                      autoAlphaRange=.1*(1:9),
                      autoAlphaDimReduceFlavor='auto',
                      summarizedExprrimentAssay=1,
                      is.counts=TRUE,
                      numcores=1,
                      ...) {
    autoAlphaMethod <- match.arg(autoAlphaMethod)

    if(class(x)=='matrix') {
        expr <- x
        se <- SummarizedExperiment::SummarizedExperiment(x)
    }
    else if(class(x)=='SummarizedExperiment') {
        expr <- SummarizedExperiment::assays(x)[[summarizedExprrimentAssay]]
        se <- x
    }
    else stop("x must be either a matrix or a SummarizedExperiment object")

    if(is.numeric(alpha)) {
        if(alpha<0 | alpha > 1) stop('alpha must be between 0 and 1')
        expr.smoothed <- smoothAndRecombine(expr, adjMatrix, alpha)
    }
    else if(alpha=='auto') {
        if(autoAlphaDimReduceFlavor=='auto') {
            autoAlphaDimReduceFlavor <- pickDimReduction(expr,
                                                         is.counts=is.counts)
            cat(paste0("Picked dimReduceFlavor: ", autoAlphaDimReduceFlavor,
                       "\n"))
        }

        smoothed.expression.matrices <- mclapply(autoAlphaRange, function(a) {
            smoothAndRecombine(expr, adjMatrix, a)
        },
        mc.cores=numcores)

        scores <- unlist(mclapply(1:length(smoothed.expression.matrices),
                                  function(i) {
            x.sm <- smoothed.expression.matrices[[i]]
            scoreSmoothing(x=x.sm, method=autoAlphaMethod,
                           is.counts=is.counts,
                           dimReduceFlavor=autoAlphaDimReduceFlavor, ...)
        },
        mc.cores=numcores))
        expr.smoothed <- smoothed.expression.matrices[[which.max(scores)]]
        chosen.a <- autoAlphaRange[which.max(scores)]
        cat(paste0("Picked alpha=",chosen.a,"\n"))
    }

    if(class(x)=='matrix') {
        return(expr.smoothed)
    }
    else if(class(x)=='SummarizedExperiment') {
        ret <- SummarizedExperiment::SummarizedExperiment(expr.smoothed,
                                                          colData=colData(x))
        return(ret)
    }
}
