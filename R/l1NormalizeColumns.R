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
