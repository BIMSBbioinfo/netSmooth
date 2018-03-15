#' Get lower dimension embedding
#'
#' @param x    gene expresison matrix [GENES x SAMPLES]
#' @param flavor    the algorithm to use to obtain the dimensionality reduction
#'                  must be in c('pca', 'tsne')
#' @param k    the number of dimensions in the reduced dimension representation
#' @param is.counts    logical: is `x` counts data
#' @param ntop    number of most variable genes to use for dimensionality
#'                reduction
#' @return    reduced dimensionality representation
#' @keywords internal
#' @importFrom scater plotPCA plotTSNE calculateCPM
#' @importFrom SingleCellExperiment reducedDim
#' @import SingleCellExperiment
dimReduce <- function(x, flavor=c('pca', 'tsne'), k=2, is.counts=TRUE, ntop=500) {
    flavor <- match.arg(flavor)
    if(flavor=='pca') function.to.call <- runPCA
    if(flavor=='tsne') function.to.call <- runTSNE

    if(is.counts){
        sce <- SingleCellExperiment(assays=list(counts=x))
        exprs(sce) <- log2(calculateCPM(sce, use_size_factors = FALSE) + 1)
    } else {
        sce <- SingleCellExperiment(assays=list(logcounts=x))
    }

    sce <- function.to.call(sce, ncomponents=k, ntop=ntop)
    red.dim <- data.frame(reducedDim(sce))[,seq_len(k)]
    colnames(red.dim) <- paste0(flavor, seq_len(k))
    return(red.dim)
}
