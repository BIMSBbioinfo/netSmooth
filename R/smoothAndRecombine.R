#' Perform network smoothing on network when the network genes and the
#' experiment genes aren't exactly the same.
#'
#' The gene network might be defined only on a subset of genes that are
#' measured in any experiment. Further, an experiment might not measure all
#' genes that are present in the network. This function projects the experiment
#' data onto the gene space defined by the network prior to smoothing. Then,
#' it projects the smoothed data back into the original dimansions.
#'
#' @param gene_expression  gene expession data to be smoothed
#'                         [N_genes x M_samples]
#' @param adj_matrix  adjacenty matrix of network to perform smoothing over.
#'                    Will be column-normalized.
#'                    Rownames and colnames should be genes.
#' @param alpha  network smoothing parameter (1 - restart probability in random
#'                walk model.
#' @param smoothing.function  must be a function that takes in data, adjacency
#'                            matrix, and alpha. Will be used to perform the
#'                            actual smoothing.
#' @usage  smoothAndRecombine(gene_expression, adj_matrix, alpha)
#' @return  matrix with network-smoothed gene expression data. Genes that are
#'          not present in smoothing network will retain original values.
#' @keywords internal
smoothAndRecombine <- function(gene_expression, adj_matrix, alpha,
                               smoothing.function=randomWalkBySolve) {
    gene_expression_in_A_space <- projectOnNetwork(gene_expression,
                                                   rownames(adj_matrix))
    gene_expression_in_A_space_smooth <- smoothing.function(
        gene_expression_in_A_space, adj_matrix, alpha)
    gene_expression_smooth <- projectFromNetworkRecombine(
        gene_expression, gene_expression_in_A_space_smooth)
    return(gene_expression_smooth)
}
