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
            fill_matrix <- Matrix::Matrix(rep(missing.value,
                prod(dim(data_in_new_space[genes_only_in_network,]))),
                nrow=length(genes_only_in_network))
        data_in_new_space[genes_only_in_network,] <- fill_matrix
        }
        return(data_in_new_space)
    }
)

setMethod("projectOnNetwork",
          signature(gene_expression='DelayedMatrix'),
          function(gene_expression, new_features, missing.value=1) {
            
            # new matrix that then will be written on disk
            data_in_new_space = as(Matrix::Matrix(rep(1, length(new_features)*
                                                     dim(gene_expression)[2]),
                                               nrow=length(new_features)), "HDF5Array")
            
            message("Type of data in new space: ", typeof(data_in_new_space@seed),"\n")
            
            rownames(data_in_new_space) <- new_features
            colnames(data_in_new_space) <- colnames(gene_expression)
            
            genes_in_both <- intersect(rownames(data_in_new_space),
                                       rownames(gene_expression))
            
            data_in_new_space[genes_in_both,] <- gene_expression[genes_in_both,]
            
            # not sure if necessary
            # genes_only_in_network <- setdiff(new_features, rownames(gene_expression))
            # if(length(genes_only_in_network)>0) {
            #   fill_matrix <- Matrix::Matrix(rep(missing.value,
            #                                     prod(dim(data_in_new_space[genes_only_in_network,]))),
            #                                 nrow=length(genes_only_in_network))
            #   data_in_new_space[genes_only_in_network,] <- fill_matrix
            # }
            return(data_in_new_space)
            
            
            
            
            
            
            
            # get genes in expression matrix and network
            in.both <- intersect(gene_expression@seed@dimnames[[1]], new_features)
            
            # get genes exclusive in network
            ppi.exclusive <- new_features[!(new_features %in% in.both)]
            
            
            expression.both <- as(gene_expression[gene_expression@seed@dimnames[[1]] %in% new_features], "HDF5Array")
            
            # init delayed matrix for ppi exclusive genes
            ppi.genes <- as( matrix(rep(missing.value,
                                        length(ppi.exclusive)*dim(gene_expression)[2]),
                                    nrow = length(ppi.exclusive)), "HDF5Array")
            
            # add col and row names
          
            #colnames(ppi.genes) <- colnames(gene_expression)
            
            # file path for ppi.genes
            #ppi.genes.filepath <- ppi.genes@seed@filepath
            
            data_in_new_space <- arbind(expression.both , ppi.genes)
            
            
            
            # add row and colnames
            colnames(data_in_new_space) <- colnames(gene_expression)
            rownames(data_in_new_space) <- c(in.both, ppi.exclusive)
            #rownames(data_in_new_space) <- 
            
            # # seed is now concat of overlap and excusive genes
            # data_in_new_space <- DelayedArray(rbind(matrix(gene_expression[
            #   gene_expression@seed@dimnames[[1]] %in% new_features],
            #   nrow = length(new_features)),
            #   matrix(rep(missing.value,
            #              length(ppi.exclusive)*dim(gene_expression)[2]),
            #          nrow = length(ppi.exclusive))
            #   ))
              
              
            # delete temporary ppi exclusive object and file on disk
            #rm(ppi.genes)
            #file.remove(ppi.genes.filepath)

            return(data_in_new_space)
          }
)
