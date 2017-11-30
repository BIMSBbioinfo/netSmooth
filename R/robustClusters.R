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
#'                           'dm' or 'auto', which uses shannon entropy to
#'                           pick the algorithm.
#' @param is.counts    logical: is the data counts
#' @param ...    arguments passed on to `clusterExperimentWorkflow`
#' @return list(clusters, proportion.robust)
#' @examples
#' x <- cbind(matrix(rexp(60000, rate=.1), ncol=100) + 1000*rexp(600, rate=.9),
#'            matrix(rexp(30000, rate=.5), ncol=50) + 1*rexp(600, rate=.9))
#' robustClusters(x)
#' @export
setMethod("robustClusters",
          signature(x='SummarizedExperiment'),
          function(x, dimReduceFlavor='auto', is.counts=TRUE, ...) {

        if(dimReduceFlavor=='auto') {
            dimReduceFlavor <- pickDimReduction(assay(x),
                                                flavors=c('pca', 'tsne'),
                                                is.counts=is.counts)
            cat(paste0("Picked dimReduceFlavor: ",dimReduceFlavor,"\n"))
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
          function(x,...) {
              robustClusters(SummarizedExperiment(x))
  }
)
