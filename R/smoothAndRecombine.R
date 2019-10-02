setGeneric(
  name = "smoothAndRecombine",
  def = function(gene_expression, adj_matrix, alpha,
                 smoothing.function=randomWalkBySolve,
                 normalizeAdjMatrix=c('rows','columns'),
                 filepath = NULL) {
    standardGeneric("smoothAndRecombine")
  }
)


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
#' @param normalizeAdjMatrix    which dimension (rows or columns) should the
#'                              adjacency matrix be normalized by. rows
#'                              corresponds to in-degree, columns to
#'                              out-degree.
#' @param filepath      String: Path to location where hdf5 output file is supposed to be saved. 
#'                      Will be ignored when regular matrices or SummarizedExperiment are
#'                      used as input.
#' @return  matrix with network-smoothed gene expression data. Genes that are
#'          not present in smoothing network will retain original values.
#' @keywords internal
setMethod("smoothAndRecombine",
          signature(gene_expression='matrix'),
          function(gene_expression, adj_matrix, alpha,
                   smoothing.function=randomWalkBySolve,
                   normalizeAdjMatrix=c('rows','columns')) {
            normalizeAdjMatrix <- match.arg(normalizeAdjMatrix)
            gene_expression_in_A_space <- projectOnNetwork(gene_expression,
                                                           rownames(adj_matrix))
            gene_expression_in_A_space_smooth <- smoothing.function(
              gene_expression_in_A_space, adj_matrix, alpha, normalizeAdjMatrix)
            gene_expression_smooth <- projectFromNetworkRecombine(
              gene_expression, gene_expression_in_A_space_smooth)
            return(gene_expression_smooth)
          })

setMethod("smoothAndRecombine",
          signature(gene_expression='Matrix'),
          function(gene_expression, adj_matrix, alpha,
                   smoothing.function=randomWalkBySolve,
                   normalizeAdjMatrix=c('rows','columns')) {
            normalizeAdjMatrix <- match.arg(normalizeAdjMatrix)
            gene_expression_in_A_space <- projectOnNetwork(gene_expression,
                                                           rownames(adj_matrix))
            gene_expression_in_A_space_smooth <- smoothing.function(
              gene_expression_in_A_space, adj_matrix, alpha, normalizeAdjMatrix)
            gene_expression_smooth <- projectFromNetworkRecombine(
              gene_expression, gene_expression_in_A_space_smooth)
            return(gene_expression_smooth)
          })

setMethod("smoothAndRecombine",
          signature(gene_expression='DelayedMatrix'),
          function(gene_expression, adj_matrix, alpha,
                   smoothing.function=randomWalkByMatrixInv,
                   normalizeAdjMatrix=c('rows','columns'),
                   filepath = NULL) {
            gene_expression_in_A_space <- projectOnNetwork(gene_expression,
                                                           rownames(adj_matrix))
            gene_expression_in_A_space_smooth <- smoothing.function(
              gene_expression_in_A_space, adj_matrix, alpha, normalizeAdjMatrix='columns')
            gene_expression_smooth <- projectFromNetworkRecombine(
              gene_expression, gene_expression_in_A_space_smooth, filepath)
            return(gene_expression_smooth)
          })


