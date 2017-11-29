#' Smooth data on graph by solving the linear equation
#' (I - alpha*A)^T * E_sm^T = E^T * (1-alpha)
#'
#'
#' @param E      initial data matrix [NxM]
#' @param A      adjacency matrix of graph to network smooth on
#'               will be column-normalized.
#' @param alpha  smoothing coefficient (1 - restart probability of
#'                random walk)
#' @return network-smoothed gene expression
#' @keywords internal
randomWalkBySolve <- function(E, A, alpha,
                              normalizeAjdMatrix=c('rows','columns')) {
    normalizeAjdMatrix <- match.arg(normalizeAjdMatrix)
    if(normalizeAjdMatrix=='rows') Anorm <- l1NormalizeRows(A)
    else if(normalizeAjdMatrix=='columns') Anorm <- l1NormalizeColumns(A)
    eye <- diag(dim(A)[1])
    AA <- t(eye - alpha*Anorm)
    BB <- (1-alpha) * E
    return(solve(AA, BB))
}
