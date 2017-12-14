#' Human Protein-Protein interaction graph
#'
#' An adjacency matrix of the 10 percent highest confidence interactions
#' between human proteins on STRINGdb
#'
#' @format A square matrix where A_{ij}=1 if gene i interacts with gene j
#' @source \url{http://www.string-db.org/}
"human.ppi"

#' Mouse Protein-Protein interaction graph
#'
#' An adjacency matrix of the 10 percent highest confidence interactions
#' between mouse proteins on STRINGdb
#'
#' @format A square matrix where A_{ij}=1 if gene i interacts with gene j
#' @source \url{http://www.string-db.org/}
"mouse.ppi"

#' A small single cell RNA-seq dataset for use in examples.
#'
#' Contains scRNAseq profiles of human blastomeres.
#'
#' @format SingleCellExperiment
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE44183}
"smallscRNAseq"

#' A small human Protein-Protein interaction graph for use in examples.
#'
#' Contains a synthetic PPI of human genes.
"smallPPI"
