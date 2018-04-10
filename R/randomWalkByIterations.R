#' Smooth data on graph by computing iterations
#'
#'
#' @param f0      initial data matrix [NxM]
#' @param adjMatrix       adjacency matrix of graph to network smooth on
#'                will be column-normalized.
#' @param alpha  smoothing coefficient (1 - restart probability of
#'                random walk)
#' @param tol    the tolerance (stopping criterion)
#' @param max.iter    the maximum number of iterations before terminating
#' @return network-smoothed gene expression
#' @keywords internal
randomWalkByIterations <- function(f0, adjMatrix, alpha,
    normalizeAjdMatrix=c('rows','columns'),
    tol=1e-6,
    max.iter=100) {
    normalizeAjdMatrix <- match.arg(normalizeAjdMatrix)
    if(normalizeAjdMatrix=='rows') Anorm <- l1NormalizeRows(adjMatrix)
    else if(normalizeAjdMatrix=='columns') Anorm <-
        l1NormalizeColumns(adjMatrix)
    f <- f0
    for(i in seq_len(max.iter)) {
        f_next <- alpha * Anorm %*% f + (1-alpha)*f0
        e <- norm(f_next - f)
        if(e<tol) break
        f <- f_next
    }
    if(i==max.iter) stop("did not converge. try increasing max.iter.")
    return(f_next)
}
