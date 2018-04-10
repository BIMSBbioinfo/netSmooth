#' Combine gene expression from smoothed space (that of the network) with the
#' expression of genes that were not smoothed (not present in network)
#' @keywords internal
#' @param original_expression    the non-smoothed expression
#' @param smoothed_expression    the smoothed gene expression, in the space
#'                               of the genes defined by the network
#' @return a matrix in the dimensions of original_expression, where values that
#'         are present in smoothed_expression are copied from there.
#' @importFrom data.table copy
projectFromNetworkRecombine <- function(original_expression,
    smoothed_expression) {
    data_in_original_space <- copy(original_expression)
    genes_in_both <- intersect(rownames(original_expression),
        rownames(smoothed_expression))
    data_in_original_space[genes_in_both,] <- as.matrix(
        smoothed_expression[genes_in_both,])
    return(data_in_original_space)
}
