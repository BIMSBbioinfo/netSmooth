clusterOne <- function(x, algorithm=c('kmeans', 'pam'), k=5) {
    algorithm <- match.arg(algorithm)
    if(algorithm == 'kmeans') yhat <- kmeans(x, k)$cluster
    else if(algorithm == 'pam') yhat <- cluster::pam(x, k)$clustering
    return(yhat)
}

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
#' @return    cluster assignments
#' @export
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
    ce <- clusterExperiment::clusterMany(se, clusterFunction=cluster.function,
                                         ks=cluster.ks,
                      isCount=is.counts, dimReduce=c("var"), nVarDims=nVarDims,
                      run=TRUE,
                      subsampling=FALSE,
                      random.seed=random.seed)

    # Run dim reduce clusterings
    if(dimReduceFlavor=='pca') dim.reduce.ks <- c(5,15,50)
    else dim.reduce.ks <- c(2,3)

    dim.reduce.clustermatrix <- matrix(rep(-1, dim(ce)[2] * length(cluster.ks) *
                                               length(dim.reduce.ks)),
                                       ncol=length(cluster.ks)*
                                           length(dim.reduce.ks))
    dim.reduce.cluster.labels <- rep(NA,
                                     length(dim.reduce.ks) * length(cluster.ks))

    i <- 1
    for(dim.reduce.k in dim.reduce.ks) {
        x <- dimReduce(assay(se), flavor=dimReduceFlavor, k=dim.reduce.k)
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
    ce <- clusterExperiment::clusterExperiment(se, cbind(clusterMatrix(ce),
                                      dim.reduce.clustermatrix),
                            transformation=transformation(ce))

    # Run a few combinations of the manys
    ce <- clusterExperiment::combineMany(ce, proportion=combineManyProportion,
                      minSize=combineManyMinSize,
                      clusterLabel="combineMany",
                      whichClusters = 'all')
    if(runMergeClusters) {
        # Make a dendrogram
        ce <- clusterExperiment::makeDendrogram(ce,dimReduce="var",ndims=500,
                             whichCluster="combineMany")

        # Merge clusters
        ce <- clusterExperiment::mergeClusters(ce, cutoff=0.05,
                                               mergeMethod="adjP",
                                               plotType="none",
                                               clusterLabel="mergeClusters")

        # Set merged to final
        ce <- clusterExperiment::setToFinal(ce,
                                            whichCluster="mergeClusters",
                                            clusterLabel="Final Clustering")
    }
    else {
        ce <- clusterExperiment::setToFinal(ce,
                                            whichCluster='combineMany',
                                            clusterLabel='Final Clustering')
    }

    # Return cluster assignments
    final.ix <- which(colnames(clusterMatrix(ce))=='Final Clustering')
    yhat <- clusterMatrix(ce)[,final.ix]
    names(yhat) <- colnames(se)
    return(yhat)
}

#' Perform robust clustering on dataset, and calculate the proportion of
#' samples in robust clusters
#' @param se    SummarizedExperiment object
#' @param dimReduceFlavor    algorithm for dimensionality reduction step
#'                           of clustering procedure. May be 'pca', 'tsne',
#'                           'dm' or 'auto', which uses shannon entropy to
#'                           pick the algorithm.
#' @param is.counts    logical: is the data counts
#' @param ...    arguments passed on to `clusterExperimentWorkflow`
#' @return list(clusters, proportion.robust)
#' @export
robustClusters <- function(se, dimReduceFlavor='auto', is.counts=TRUE, ...) {
    if(dimReduceFlavor=='auto') {
        dimReduceFlavor <- pickDimReduction(assay(se),
                                flavors=c('pca', 'tsne'), is.counts=is.counts)
        cat(paste0("Picked dimReduceFlavor: ",dimReduceFlavor,"\n"))
    }
    yhat <- clusterExperimentWorkflow(se, is.counts=is.counts,
                                      dimReduceFlavor=dimReduceFlavor, ...)
    proportion.robust <- mean(yhat!=-1)
    return(list(clusters=yhat, proportion.robust=proportion.robust))
}
