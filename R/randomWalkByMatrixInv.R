setGeneric(
  name = "randomWalkByMatrixInv",
  def = function(f0, adjMatrix, alpha,
                 normalizeAdjMatrix=c('rows','columns')) {
    standardGeneric("randomWalkByMatrixInv")
  }
)

#' Smooth data on graph by computing the closed-form steady state
#' distribution of the random walk with restarts process.
#'
#' The closed-form solution is given by
#'   f_{ss} = (1 - alpha) * (I - alpha * A)^{-1} * f_0
#' and is computed by matrix inversion in this function.
#'
#' @param f0      initial data matrix [NxM]
#' @param adjMatrix       adjacency matrix of graph to network smooth on
#'                will be column-normalized.
#' @param alpha  smoothing coefficient (1 - restart probability of
#'                random walk)
#' @return network-smoothed gene expression
#' @keywords internal
setMethod("randomWalkByMatrixInv",
          signature(f0='matrix'),
          function(f0, adjMatrix, alpha,
          normalizeAdjMatrix=c('rows','columns')) {
            normalizeAdjMatrix <- match.arg(normalizeAdjMatrix)
            if(normalizeAdjMatrix=='rows') Anorm <- l1NormalizeRows(adjMatrix)
            else if(normalizeAdjMatrix=='columns') Anorm <-
                l1NormalizeColumns(adjMatrix)
            eye <- diag(dim(adjMatrix)[1])
            K <- (1 - alpha) * solve(eye - alpha * Anorm)
            return(K %*% f0)
            })

setMethod("randomWalkByMatrixInv",
          signature(f0='Matrix'),
          function(f0, adjMatrix, alpha,
                   normalizeAdjMatrix=c('rows','columns')) {
            normalizeAdjMatrix <- match.arg(normalizeAdjMatrix)
            if(normalizeAdjMatrix=='rows') Anorm <- l1NormalizeRows(adjMatrix)
            else if(normalizeAdjMatrix=='columns') Anorm <-
                l1NormalizeColumns(adjMatrix)
            eye <- diag(dim(adjMatrix)[1])
            K <- (1 - alpha) * solve(eye - alpha * Anorm)
            return(K %*% f0)
          })

setMethod("randomWalkByMatrixInv",
          signature(f0='DelayedMatrix'),
          function(f0, adjMatrix, alpha,
                   normalizeAdjMatrix) {
            normalizeAdjMatrix <- match.arg(normalizeAdjMatrix)
            if(normalizeAdjMatrix=='rows') Anorm <- l1NormalizeRows(adjMatrix)
            else if(normalizeAdjMatrix=='columns') Anorm <-
                l1NormalizeColumns(adjMatrix)
            eye <- diag(dim(adjMatrix)[1])
            K <- (1 - alpha) * solve(eye - alpha * Anorm)
            return(K %*% f0)
          })

