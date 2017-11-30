setGeneric(
    name = "pickDimReduction",
    def = function(x, ...) {
        standardGeneric("pickDimReduction")
    }
)

#' Pick the dimensionality reduction method for a dataset that gives the
#' 2D embedding with the highest entropy
#'
#' @rdname pickDimReduction
#' @aliases pickDimReduction
#' @param x  matrix or SummarizedExperiment object [GENES x SAMPLES]
#' @param flavors    list of dimensionality reduction algorithms to try.
#'                   Currently the options are "pca" and "tsne"
#' @param is.counts    logical: is exprs count data
#' @return    name of dimensionality reduction method that gives the highest
#'            2d entropy
#' @importFrom SummarizedExperiment assay
#' @examples
#' x <- matrix(rnbinom(60000, size=1, prob = .1), ncol=100)
#' pickDimReduction(x)
#' @export
setMethod("pickDimReduction",
          signature(x='matrix'),
          function(x, flavors=c('pca', 'tsne'), is.counts=TRUE) {
        if(class(x)=='matrix') x <- x
        else if(class(x)=='SummarizedExperiment') x <- assay(x)
        else stop("must be matrix or SummarizedExperiment")
        entropies <- sapply(flavors, function(flavor) {
            calc2DEntropy(dimReduce(x, flavor=flavor, is.counts=is.counts))
        })
        return(names(which.max(entropies)))
    }
)

#' @export
#' @importFrom SummarizedExperiment assay
#' @rdname pickDimReduction
setMethod("pickDimReduction",
          signature(x='SummarizedExperiment'),
          function(x) pickDimReduction(assay(x)))
