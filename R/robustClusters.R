setGeneric(
    name = "robustClusters",
    def = function(x, ...) {
        standardGeneric("robustClusters")
    }
)

#' Perform robust clustering on dataset, and calculate the proportion of
#' samples in robust clusters
#' @rdname robustClusters
#' @aliases robustClusters
#' @param x    matrix or SummarizedExperiment object
#' @param dimReduceFlavor    algorithm for dimensionality reduction step
#'                           of clustering procedure. May be 'pca', 'tsne',
#'                           'dm', 'umap' or 'auto', which uses shannon entropy to
#'                           pick the algorithm.
#' @param is.counts    logical: is the data counts
#' @param ...    arguments passed on to `clusterExperimentWorkflow`
#' @examples
#' data("smallscRNAseq")
#' robustClusters(smallscRNAseq, dimReduceFlavor='pca')
#' @return list(clusters, proportion.robust)
#' @export
setMethod("robustClusters",
    signature(x='SummarizedExperiment'),
    function(x, dimReduceFlavor='auto', is.counts=TRUE, ...) {
        if(any(is(assay(x))=='Matrix')){
            stop("robust clustering does not currently support sparse matrices.")
        }
        if(dimReduceFlavor=='auto') {
            dimReduceFlavor <- pickDimReduction(assay(x),
                flavors=c('pca', 'tsne', 'umap'),
                is.counts=is.counts)
            message("Picked dimReduceFlavor: ",dimReduceFlavor,"\n")
        }
        yhat <- clusterExperimentWorkflow(x, is.counts=is.counts,
            dimReduceFlavor=dimReduceFlavor, ...)
        proportion.robust <- mean(yhat!=-1)
        return(list(clusters=yhat, proportion.robust=proportion.robust))
    }
)

#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @rdname robustClusters
setMethod("robustClusters",
    signature(x='matrix'),
    function(x, ...) {
        robustClusters(SummarizedExperiment(x), ...)
    }
)

