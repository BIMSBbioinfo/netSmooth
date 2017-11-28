projectFromNetworkRecombine <- function(original_expression,
                                        smoothed_expression) {
    data_in_original_space <- copy(original_expression)
    genes_in_both <- intersect(rownames(original_expression),
                               rownames(smoothed_expression))
    data_in_original_space[genes_in_both,] <- as.matrix(
        smoothed_expression[genes_in_both,])
    return(data_in_original_space)
}
