#' Smooth data on graph by computing the closed-form steady state
#' distribution of the random walk with restarts process.
#'
#' The closed-form solution is given by
#'   f_{ss} = (1 - alpha) * (I - alpha * A)^{-1} * f_0
#' and is computed by matrix inversion in this function.
#'
#' @param f0      initial data matrix [NxM]
#' @param adj_matrix       adjacency matrix of graph to network smooth on
#'                will be column-normalized.
#' @param alpha  smoothing coefficient (1 - restart probability of
#'                random walk)
#' @usage randomWalkByMatrixInv(gene_expression, adj_matrix, alpha)
#' @return network-smoothed gene expression
#' @keywords internal
randomWalkByMatrixInv <- function(f0, A, alpha,
                                  normalizeAjdMatrix=c('rows','columns')) {
    normalizeAjdMatrix <- match.arg(normalizeAjdMatrix)
    if(normalizeAjdMatrix=='rows') Anorm <- l1NormalizeRows(A)
    else if(normalizeAjdMatrix=='columns') Anorm <- l1NormalizeColumns(A)
    eye <- diag(dim(A)[1])
    K <- (1 - alpha) * solve(eye - alpha * Anorm)
    return(K %*% f0)
}
