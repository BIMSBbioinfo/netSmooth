#' Column-normalize a sparse, symmetric matrix (using the l1 norm) so that each
#' column sums to 1.
#'
#' @param A matrix
#' @usage l1NormalizeColumns(A)
#' @return row-normalized sparse matrix object
#' @keywords internal
l1NormalizeColumns <- function(A) {
    return(Matrix::t(Matrix::t(A)/Matrix::colSums(A)))
}

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
randomWalkByMatrixInv <- function(f0, A, alpha) {
    Anorm <- l1NormalizeColumns(A)
    eye <- diag(dim(A)[1])
    K <- (1 - alpha) * solve(eye - alpha * Anorm)
    return(K %*% f0)
}

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
randomWalkBySolve <- function(E, A, alpha) {
    E <- t(E)
    Anorm <- l1NormalizeColumns(A)
    eye <- diag(dim(A)[1])
    AA <- t(eye - alpha*Anorm)
    BB <- (1-alpha) * t(E)
    return(solve(AA, BB))
}
