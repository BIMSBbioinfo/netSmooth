setGeneric(
    name = "projectOnNetwork",
    def = function(gene_expression, new_features, missing.value=0) {
        standardGeneric("projectOnNetwork")
    }
)

#' Project the gene expression matrix onto a lower space
#' of the genes defined in the smoothing network
#' @param gene_expression    gene expression matrix
#' @param new_features       the genes in the network, on which to project
#'                           the gene expression matrix
#' @param missing.value      value to assign to genes that are in network,
#'                           but missing from gene expression matrix
#' @return the gene expression matrix projected onto the gene space defined by new_features
#' @rdname projectOnNetwork
#' @aliases projectOnNetwork
#' @keywords internal
setMethod("projectOnNetwork",
    signature(gene_expression='matrix'),
    function(gene_expression, new_features, missing.value=0) {
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
})

setMethod("projectOnNetwork",
          signature(gene_expression='Matrix'),
          function(gene_expression, new_features, missing.value=0) {
              data_in_new_space = Matrix::Matrix(rep(0, length(new_features)*
                                                 dim(gene_expression)[2]),
                                         nrow=length(new_features))
              rownames(data_in_new_space) <- new_features
              colnames(data_in_new_space) <- colnames(gene_expression)
              genes_in_both <- intersect(rownames(data_in_new_space),
                                         rownames(gene_expression))
              data_in_new_space[genes_in_both,] <- gene_expression[genes_in_both,]

              genes_only_in_network <- setdiff(new_features, rownames(gene_expression))
              if(length(genes_only_in_network)>0) {
                  fill_matrix <- Matrix::Matrix(rep(missing.value, prod(dim(data_in_new_space[genes_only_in_network,]))),
                                        nrow=length(genes_only_in_network))
                  data_in_new_space[genes_only_in_network,] <- fill_matrix
              }
              return(data_in_new_space)
          })
