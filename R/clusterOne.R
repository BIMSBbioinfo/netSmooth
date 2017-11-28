clusterOne <- function(x, algorithm=c('kmeans', 'pam'), k=5) {
    algorithm <- match.arg(algorithm)
    if(algorithm == 'kmeans') yhat <- kmeans(x, k)$cluster
    else if(algorithm == 'pam') yhat <- cluster::pam(x, k)$clustering
    return(yhat)
}
