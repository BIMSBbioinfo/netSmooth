setGeneric(
  name = "projectFromNetworkRecombine",
  def = function(original_expression, smoothed_expression) {
    standardGeneric("projectFromNetworkRecombine")
  }
)

#' Combine gene expression from smoothed space (that of the network) with the
#' expression of genes that were not smoothed (not present in network)
#' @keywords internal
#' @param original_expression    the non-smoothed expression
#' @param smoothed_expression    the smoothed gene expression, in the space
#'                               of the genes defined by the network
#' @return a matrix in the dimensions of original_expression, where values that
#'         are present in smoothed_expression are copied from there.
#' @importFrom data.table copy
setMethod("projectFromNetworkRecombine",
          signature(original_expression='matrix'),
          function(original_expression,
                   smoothed_expression) {
            data_in_original_space <- copy(original_expression)
            genes_in_both <- intersect(rownames(original_expression),
                                       rownames(smoothed_expression))
            data_in_original_space[genes_in_both,] <- as.matrix(
              smoothed_expression[genes_in_both,])
            return(data_in_original_space)
          })

setMethod("projectFromNetworkRecombine",
          signature(original_expression='Matrix'),
          function(original_expression,
                   smoothed_expression) {
            data_in_original_space <- copy(original_expression)
            genes_in_both <- intersect(rownames(original_expression),
                                       rownames(smoothed_expression))
            data_in_original_space[genes_in_both,] <- as.matrix(
              smoothed_expression[genes_in_both,])
            return(data_in_original_space)
          })

setMethod("projectFromNetworkRecombine",
          signature(original_expression='DelayedMatrix'),
          function(original_expression,
                   smoothed_expression) {
            
            
            # genes in orignal and smoothed expression
            genes.in.both <- intersect(original_expression@seed@dimnames[[1]],
                      smoothed_expression@seed@dimnames[[1]])
            
            # genes not in intersection
            original.exclusive <- original_expression@seed@dimnames[[1]][which(!(original_expression@seed@dimnames[[1]] %in% genes.in.both))]
              

            # new File, empty matrix as seed
            data_in_original_space <- as( matrix(rep(0, dim(original_expression)[1]*dim(original_expression)[2]),
                                                 nrow=dim(original_expression)[1]), "HDF5Array")
            
            # set and row/col names and coerce to DelayedMatrix
            rownames(data_in_original_space) <- original_expression@seed@dimnames[[1]]
            colnames(data_in_original_space) <- original_expression@seed@dimnames[[2]]
            
            # replace rows in output
            data_in_original_space[genes.in.both,] <- smoothed_expression[genes.in.both,]
            data_in_original_space[original.exclusive,] <- original_expression[original.exclusive,]
            
            # delete file containing smoothing matrix
            file.path <- path(smoothed_expression)
            file.remove(file.path)

            return(data_in_original_space)
          })
  

