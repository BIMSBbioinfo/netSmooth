projectOnNetwork <- function(gene_expression, new_features, missing.value=0) {
    data_in_new_space = matrix(rep(0, length(new_features)*
                                       dim(gene_expression)[2]),
                               nrow=length(new_features))
    rownames(data_in_new_space) <- new_features
    colnames(data_in_new_space) <- colnames(gene_expression)
    genes_in_both <- intersect(rownames(data_in_new_space),
                               rownames(gene_expression))
    data_in_new_space[genes_in_both,] <- gene_expression[genes_in_both,]

    genes_only_in_network <- setdiff(new_features, rownames(gene_expression))
    data_in_new_space[genes_only_in_network,] <- missing.value
    return(data_in_new_space)
}
