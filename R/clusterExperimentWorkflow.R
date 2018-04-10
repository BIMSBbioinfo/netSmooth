#' Performs clustering workflow using `clusterExperiment` functions
#' @param se    SummarizedExperiment object
#' @param dimReduceFlavor    algorithm for reduced dimension embedding step
#' @param cluster.ks    range of Ks to cluster over
#' @param cluster.function    clustering algorithm to use for all clusterings
#' @param nVarDims    numbers of variable genes to perform clusterings over
#' @param combineManyProportion    proportion of times samples need to be
#'                                 co-clustered for co-clustering step
#' @param combineManyMinSize    minimum cluster size
#' @param runMergeClusters    logical: merge similar clusters
#' @param is.counts    logical: is data counts
#' @param random.seed    passed to clusterExperiment. set to NULL in order
#'                       to not set a random seed.
#' @return    cluster assignments
#' @importFrom SummarizedExperiment assay
#' @importFrom clusterExperiment clusterMany clusterMatrix ClusterExperiment
#'                               combineMany makeDendrogram mergeClusters
#'                               setToFinal transformation
#' @keywords internal
clusterExperimentWorkflow <- function(se,
    dimReduceFlavor=c('pca', 'tsne', 'dm'),
    cluster.ks=5:10,
    cluster.function='pam',
    nVarDims=c(100,500,1000),
    combineManyProportion=.7,
    combineManyMinSize=4,
    runMergeClusters=TRUE,
    is.counts=TRUE,
    random.seed=1) {
    dimReduceFlavor <- match.arg(dimReduceFlavor)
    # Run variable genes clusterings
    nVarDims <- nVarDims[nVarDims < dim(se)[1]]
    ce <- clusterMany(se, clusterFunction=cluster.function,
        ks=cluster.ks,
        isCount=is.counts, dimReduce=c("var"), nVarDims=nVarDims,
        run=TRUE,
        subsampling=FALSE,
        random.seed=random.seed)

    # Run dim reduce clusterings
    if(dimReduceFlavor=='pca') dim.reduce.ks <- c(5,15,50)
    else dim.reduce.ks <- c(2,3)

    dim.reduce.ks <- dim.reduce.ks[dim.reduce.ks <= dim(se)[1]]
    dim.reduce.ks <- dim.reduce.ks[dim.reduce.ks <= dim(se)[2]]

    dim.reduce.clustermatrix <- matrix(rep(-1, dim(ce)[2] * length(cluster.ks) *
        length(dim.reduce.ks)),
        ncol=length(cluster.ks)*
        length(dim.reduce.ks))
    dim.reduce.cluster.labels <- rep(NA,
        length(dim.reduce.ks) * length(cluster.ks))

    i <- 1
    for(dim.reduce.k in dim.reduce.ks) {
        x <- dimReduce(assay(se), flavor=dimReduceFlavor, k=dim.reduce.k,
            is.counts=is.counts)
        for(cluster.k in cluster.ks) {
            dim.reduce.cluster.labels[i] <- paste0(dimReduceFlavor,
                dim.reduce.k,
                cluster.function,
                cluster.k)
            dim.reduce.clustermatrix[,i] <- clusterOne(x, cluster.function,
                cluster.k)
            i <- i+1
        }
    }
    colnames(dim.reduce.clustermatrix) <- dim.reduce.cluster.labels

    # Make a new overall clusterExperiment object
    ce <- ClusterExperiment(se,
        cbind(clusterMatrix(ce),
        dim.reduce.clustermatrix),
        transformation=transformation(ce))

    # Run a few combinations of the manys
    ce <- combineMany(ce, proportion=combineManyProportion,
        minSize=combineManyMinSize,
        clusterLabel="combineMany",
        whichClusters = 'all')
    if(runMergeClusters) {
        # Make a dendrogram
        ce <- makeDendrogram(ce, reduceMethod="var", nDims=500,
            whichCluster="combineMany")

        # Merge clusters
        ce <- mergeClusters(ce, cutoff=0.05,
            mergeMethod="adjP",
            plotInfo="none",
            plot=FALSE,
            showWarnings=FALSE,
            clusterLabel="mergeClusters")

        # Set merged to final
        ce <- setToFinal(ce,
            whichCluster="mergeClusters",
            clusterLabel="Final Clustering")
    }
    else {
        ce <- setToFinal(ce,
            whichCluster='combineMany',
            clusterLabel='Final Clustering')
    }

    # Return cluster assignments
    final.ix <- which(colnames(clusterMatrix(ce))== 'Final Clustering')
    yhat <- clusterMatrix(ce)[,final.ix]
    names(yhat) <- colnames(se)
    return(yhat)
}
