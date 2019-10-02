require(STRINGdb)
require(igraph)
require(biomaRt)

getSTRINGdbForSpecies <- function(species=c('human','mouse')) {
    species <- match.arg(species)
    if(species=='human') return( STRINGdb$new(species=9606) )
    if(species=='mouse') return( STRINGdb$new(species=10090) )

}

getBiomartForspecies <- function(species=c('human','mouse')) {
    species <- match.arg(species)
    if(species=='human') return( useMart(host = 'grch37.ensembl.org',
                                         biomart='ENSEMBL_MART_ENSEMBL',
                                         dataset='hsapiens_gene_ensembl') )
    if(species=='mouse') return( useMart("ensembl",
                                         dataset="mmusculus_gene_ensembl") )
}

getPPIFromStringDB <- function(species) {
    string_db <- getSTRINGdbForSpecies(species)
    human_graph <- string_db$get_graph()

    edge.scores <- E(human_graph)$combined_score
    ninetyth.percentile <- quantile(edge.scores, 0.9)
    thresh <- data.frame(name='90th percentile', val=ninetyth.percentile)



    human_graph <- subgraph.edges(human_graph,
                                  E(human_graph)[combined_score >
                                                     ninetyth.percentile])

    adj_matrix <- as_adjacency_matrix(human_graph)

    protein_ids <- sapply(strsplit(rownames(adj_matrix), '\\.'),
                          function(x) x[2])

    mart = getBiomartForspecies(species)

    mart_results <- getBM(attributes = c("ensembl_gene_id",
                                         "ensembl_peptide_id"),
    filters = "ensembl_peptide_id", values = protein_ids,
    mart = mart)

    ix <- match(protein_ids, mart_results$ensembl_peptide_id)
    ix <- ix[!is.na(ix)]



    newnames <- protein_ids
    newnames[match(mart_results[ix,'ensembl_peptide_id'], newnames)] <-
        mart_results[ix, 'ensembl_gene_id']
    rownames(adj_matrix) <- newnames
    colnames(adj_matrix) <- newnames

    ppi <- adj_matrix[!duplicated(newnames), !duplicated(newnames)]
    nullrows <- Matrix::rowSums(ppi)==0
    ppi <- ppi[!nullrows,!nullrows]
    return(ppi)
}

human.ppi <- getPPIFromStringDB('human')
mouse.ppi <- getPPIFromStringDB('mouse')

devtools::use_data(human.ppi, overwrite = TRUE)
devtools::use_data(mouse.ppi, overwrite = TRUE)
