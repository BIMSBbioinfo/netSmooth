#' Row-normalize a sparse, symmetric matrix (using the l1 norm) so that each
#' row sums to 1.
#'
#' @param A matrix
#' @usage l1NormalizeRows(A)
#' @return row-normalized sparse matrix object
#' @keywords internal
l1NormalizeRows <- function(A) {
    return(A/Matrix::rowSums(A))
}
