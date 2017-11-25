#' Row-normalize a sparse matrix (using the l1 norm) so that each
#' row sums to 1.
#'
#' @param A matrix
#' @usage l1NormalizeRows(A)
#' @return row-normalized sparse matrix object
l1NormalizeRows <- function(A) {
    rs = Matrix::rowSums(A)
    factors = 1 / replace(rs, rs==0, 1)
    return(Matrix::Diagonal(x=factors) %*% A)
}

#' Column-normalize a sparse, symmetric matrix (using the l1 norm) so that each
#' column sums to 1.
#'
#' @param A matrix
#' @usage l1NormalizeColumns(A)
#' @return row-normalized sparse matrix object
l1NormalizeColumns <- function(A) {
    return(Matrix::t(l1NormalizeRows((A))))
}

#' Smooth data on graph by computing the closed-form steady state
#' distribution of the random walk with restarts process.
#'
#' The closed-form solution is given by
#'   f_{ss} = (1 - lambda) * (I - lambda * A)^{-1} * f_0
#' and is computed by matrix inversion in this function.
#'
#' @param f0      initial data matrix [NxM]
#' @param adj_matrix       adjacency matrix of graph to network smooth on
#'                will be column-normalized.
#' @param lambda  smoothing coefficient (1 - restart probability of
#'                random walk)
#' @usage randomWalkByMatrixInv(gene_expression, adj_matrix, lambda)
#' @return network-smoothed gene expression
#' @export
randomWalkByMatrixInv <- function(f0, A, lambda) {
    Anorm <- l1NormalizeColumns(A)
    eye <- diag(dim(A)[1])
    K <- (1 - lambda) * solve(eye - lambda * Anorm)
    return(K %*% f0)
}
