#' Smooth data on graph by performing random walk iterations until convergence
#'
#' @param f0      initial data matrix [NxM]
#' @param adj_matrix       adjacency matrix of graph to network smooth on
#'                will be column-normalized.
#' @param lambda  smoothing coefficient (1 - restart probability of
#'                random walk)
#' @param max.iter  maximum number of iterations to allow for convergence.
#'                  if not converfed in max.iter tries, gives warning.
#' @param tol     error tolerance for convergence criterion:
#'                will terminate when norm(f_t - f_{t+1}) < tol
#' @usage random.walk.by.iterations(gene_expression, adj_matrix, lambda)
#' @return network-smoothed gene expression
#' @export
random.walk.by.iterations <- function(f0, A, lambda, max.iter=1000, tol=1e-3) {
    Anorm <- l1.normalize.columns(A)
    e <- Inf
    f <- f0
    for(k in 1:max.iter){
        f_next <- lambda * Anorm %*% f + (1-lambda) * f0
        e <- Matrix::norm(f_next-f)
        f <- f_next
        if(e < tol){
            return(f)
        }
    }
    warning("Did not converge. Make sure `A` is row-normalized. Try increasing `max.iter`.")
    return(f)
}

#' Row-normalize a sparse matrix (using the l1 norm) so that each
#' row sums to 1.
#'
#' @param A matrix
#' @usage l1.normalize.rows(A)
#' @return row-normalized sparse matrix object
l1.normalize.rows <- function(A) {
    rs = Matrix::rowSums(A)
    factors = 1 / replace(rs, rs==0, 1)
    return(Matrix::Diagonal(x=factors) %*% A)
}

#' Column-normalize a sparse, symmetric matrix (using the l1 norm) so that each
#' column sums to 1.
#'
#' @param A matrix
#' @usage l1.normalize.columns(A)
#' @return row-normalized sparse matrix object
l1.normalize.columns <- function(A) {
    return(Matrix::t(l1.normalize.rows((A))))
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
#' @usage random.walk.by.matrix.inv(gene_expression, adj_matrix, lambda)
#' @return network-smoothed gene expression
#' @export
random.walk.by.matrix.inv <- function(f0, A, lambda) {
    Anorm <- l1.normalize.columns(A)
    eye <- diag(dim(A)[1])
    K <- (1 - lambda) * solve(eye - lambda * Anorm)
    return(K %*% f0)
}
