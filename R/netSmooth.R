setGeneric(
    name = "netSmooth",
    def = function(x, ...) {
        standardGeneric("netSmooth")
    }
)

#' Perform network smoothing of gene expression or other omics data
#' @param x    matrix or SummarizedExperiment
#' @param adjMatrix    adjacency matrix of gene network to use
#' @param alpha    numeric in [0,1] or 'audo'. if 'auto', the optimal
#'                 value for alpha will be automatically chosen among the values
#'                 specified in `autoAlphaRange`, using the strategy
#'                 specified in `autoAlphaMethod`
#' @param normalizeAdjMatrix    how to normalize the adjacency matrix
#'                              possible values are 'rows' (in-degree)
#'                              and 'columns' (out-degree)
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
#' @param is.counts    logical: is the assay count data
#' @param bpparam    instance of SnowParam, for parallel computation with the
#'                   `alpha='auto'` option. See the BiocParallel manual.
#' @param ...    arguments passed on to `robustClusters` if using the robustness
#'               criterion for optimizing alpha
#' @return network-smoothed gene expression matrix or SummarizedExperiment
#'         object
#' @examples
#' x <- matrix(rnbinom(12000, size=1, prob = .1), ncol=60)
#' rownames(x) <- paste0('gene', seq_len(dim(x)[1]))
#'
#' adj_matrix <- matrix(as.numeric(rnorm(200*200)>.8), ncol=200)
#' rownames(adj_matrix) <- colnames(adj_matrix) <- paste0('gene', seq_len(dim(x)[1]))
#' x.smoothed <- netSmooth(x, adj_matrix, alpha=0.5)
#' @export
#' @rdname netSmooth
#' @inheritParams netSmooth,matrix-method
#' @aliases netSmooth
#' @importFrom SummarizedExperiment colData
setMethod("netSmooth",
    signature(x='matrix'),
    function(x, adjMatrix, alpha='auto',
        normalizeAdjMatrix=c('rows','columns'),
        autoAlphaMethod=c('robustness', 'entropy'),
        autoAlphaRange=.1*(seq_len(9)),
        autoAlphaDimReduceFlavor='auto',
        is.counts=TRUE,
        bpparam=BiocParallel::SnowParam(workers = 1, type = "SOCK"),
        ...) {
        autoAlphaMethod <- match.arg(autoAlphaMethod)
        normalizeAdjMatrix <- match.arg(normalizeAdjMatrix)

        stopifnot(is(adjMatrix, 'matrix') || is(adjMatrix, 'sparseMatrix'))
        stopifnot((is.numeric(alpha) && (alpha > 0 && alpha < 1)) || alpha == "auto")

        if(is.numeric(alpha)) {
            message("Using given alpha: ", alpha,"\n")
            if(alpha<0 | alpha > 1) {
                stop('alpha must be between 0 and 1')
            }
            x.smoothed <- smoothAndRecombine(x, adjMatrix, alpha,
                                        normalizeAdjMatrix=normalizeAdjMatrix)
        } else if(alpha=='auto') {
            if(autoAlphaDimReduceFlavor=='auto') {
                autoAlphaDimReduceFlavor <-
                    pickDimReduction(x, is.counts=is.counts)
                message("Picked dimReduceFlavor: ",
                    autoAlphaDimReduceFlavor, "\n")
            }

            smoothed.expression.matrices <- BiocParallel::bplapply(
                autoAlphaRange,
                function(a) {
                    smoothAndRecombine(x, adjMatrix, a,
                        normalizeAdjMatrix=normalizeAdjMatrix)
                }
            )

            scores <- unlist(BiocParallel::bplapply(
                seq_len(length(smoothed.expression.matrices)),
                function(i) {
                    x.sm <- smoothed.expression.matrices[[i]]
                    scoreSmoothing(x=x.sm,
                        method=autoAlphaMethod,
                        is.counts=is.counts,
                        dimReduceFlavor=autoAlphaDimReduceFlavor, ...)
                }))
            x.smoothed <- smoothed.expression.matrices[[which.max(scores)]]
            chosen.a <- autoAlphaRange[which.max(scores)]
            message("Picked alpha=",chosen.a,"\n")
        } else stop("unsupprted alpha value: ", class(alpha))
        return(x.smoothed)
    }
)

#' @rdname netSmooth
#' @export
setMethod("netSmooth",
    signature(x='SummarizedExperiment'),
    function(x, ...) {
        matrixdata <- assay(x)
        ret <- netSmooth(matrixdata, ...)
        return(SummarizedExperiment(ret, colData=colData(x)))
    })

#' @rdname netSmooth
#' @export
setMethod("netSmooth",
    signature(x='Matrix'),
    function(x, adjMatrix, alpha='auto',
        normalizeAdjMatrix=c('rows','columns'),
        autoAlphaMethod=c('robustness', 'entropy'),
        autoAlphaRange=.1*(seq_len(9)),
        autoAlphaDimReduceFlavor='auto',
        is.counts=TRUE,
        bpparam=BiocParallel::SnowParam(workers = 1, type = "SOCK"),
        ...) {
        autoAlphaMethod <- match.arg(autoAlphaMethod)
        normalizeAdjMatrix <- match.arg(normalizeAdjMatrix)

        stopifnot(is(adjMatrix, 'matrix') || is(adjMatrix, 'sparseMatrix'))
        stopifnot((is.numeric(alpha) && (alpha > 0 && alpha < 1)) || alpha == "auto")

        if(is.numeric(alpha)) {
            message("Using given alpha: ", alpha,"\n")
            if(alpha<0 | alpha > 1) {
                stop('alpha must be between 0 and 1')
            }
            x.smoothed <- smoothAndRecombine(x, adjMatrix, alpha,
            normalizeAdjMatrix=normalizeAdjMatrix)
        } else if(alpha=='auto') {
            if(autoAlphaDimReduceFlavor=='auto') {
                autoAlphaDimReduceFlavor <- pickDimReduction(x,
                    is.counts=is.counts)
                message("Picked dimReduceFlavor: ", autoAlphaDimReduceFlavor,
                    "\n")
            }

            smoothed.expression.matrices <- BiocParallel::bplapply(
                autoAlphaRange,
                function(a) {
                    smoothAndRecombine(x, adjMatrix, a,
                        normalizeAdjMatrix=normalizeAdjMatrix)
                }
            )

            scores <- unlist(BiocParallel::bplapply(
                seq_len(length(smoothed.expression.matrices)),
                function(i) {
                    x.sm <- smoothed.expression.matrices[[i]]
                    scoreSmoothing(x=x.sm,
                        method=autoAlphaMethod,
                        is.counts=is.counts,
                        dimReduceFlavor=autoAlphaDimReduceFlavor, ...)
                }
            ))
            x.smoothed <- smoothed.expression.matrices[[which.max(scores)]]
            chosen.a <- autoAlphaRange[which.max(scores)]
            message("Picked alpha=",chosen.a,"\n")
        } else stop("unsupprted alpha value: ", class(alpha))
        return(x.smoothed)
    }
)
