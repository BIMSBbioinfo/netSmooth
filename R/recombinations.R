project.on.network <- function(gene_expression, new_features, missing.value=0) {
    data_in_new_space = matrix(rep(0, length(new_features)*dim(gene_expression)[2]), nrow=length(new_features))
    rownames(data_in_new_space) <- new_features
    colnames(data_in_new_space) <- colnames(gene_expression)
    genes_in_both <- intersect(rownames(data_in_new_space), rownames(gene_expression))
    data_in_new_space[genes_in_both,] <- gene_expression[genes_in_both,]
    return(data_in_new_space)
}

project.from.network.recombine <- function(original_expression, smoothed_expression) {
    data_in_original_space <- copy(original_expression)
    genes_in_both <- intersect(rownames(original_expression), rownames(smoothed_expression))
    data_in_original_space[genes_in_both,] <- as.matrix(smoothed_expression[genes_in_both,])
    return(data_in_original_space)
}

smooth.and.recombine <- function(gene_expression, adjnorm, lambda, smoothing.function=netsmooth::random.walk) {
    gene_expression_in_A_space <- project.on.network(gene_expression, rownames(adjnorm))
    gene_expression_in_A_space_smooth <- smoothing.function(gene_expression_in_A_space, adjnorm, lambda)
    gene_expression_smooth <- project.from.network.recombine(gene_expression, gene_expression_in_A_space_smooth)
    return(gene_expression_smooth)
}
