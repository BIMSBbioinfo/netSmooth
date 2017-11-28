#' Get lower dimension embedding
#'
#' @param expr    gene expresison matrix [GENES x SAMPLES]
#' @param flavor    the algorithm to use to obtain the dimensionality reduction
#'                  must be in c('pca', 'tsne')
#' @param k    the number of dimensions in the reduced dimension representation
#' @param is.counts    logical: is `expr` counts data
#' @return    reduced dimensionality representation
#' @keywords internal
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
