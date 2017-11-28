#' Calculate entropy in 2D data
#' @param x    the 2D data to get entropy from
#' @param numBins1    the number of bins along the first dimension to
#'                    discretize data into
#' @param numBins2    the number of bins along the second dimension to
#'                    discretize data into
#' @return The Shannon entropy in the 2D data x
#' @keywords internal
calc2DEntropy <- function(x, numBins1=20, numBins2=20) {
    embedding.discretized <- entropy::discretize2d(x[,1], x[,2], numBins1, numBins2)
    return(entropy::entropy(embedding.discretized))
}
