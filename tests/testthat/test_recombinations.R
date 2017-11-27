context("Recombinations")

test_that("projectOnNetwork can subset", {
    expr <- matrix(rep(1,20), ncol=5)
    rownames(expr) <- paste('gene', 1:(dim(expr)[1]))
    network.genes <- c('gene 1', 'gene 2')
    expect_equal(projectOnNetwork(expr, network.genes), expr[network.genes,])
})

test_that("projectOnNetwork can superset", {
    expr <- matrix(rep(1,20), ncol=5)
    rownames(expr) <- paste('gene', 1:(dim(expr)[1]))
    network.genes <- c('gene 1', 'gene 5')
    projected <- projectOnNetwork(expr, network.genes)
    expect_equal(rownames(projected), network.genes)
})

test_that("projectOnNetwork fills missing values", {
    expr <- matrix(rep(1,20), ncol=5)
    rownames(expr) <- paste('gene', 1:(dim(expr)[1]))
    network.genes <- c('gene 1', 'gene 5')
    projected <- projectOnNetwork(expr, network.genes)
    expect_equal(projected['gene 5',], rep(0, dim(expr)[2]))

    projected <- projectOnNetwork(expr, network.genes, missing.value = -1)
    expect_equal(projected['gene 5',], rep(-1, dim(expr)[2]))

    projected <- projectOnNetwork(expr, network.genes, missing.value = NA)
    expect_equal(projected['gene 5',], as.numeric(rep(NA, dim(expr)[2])))
})

############################

expr <- matrix(rep(1,20), ncol=5)
rownames(expr) <- paste('gene', 1:(dim(expr)[1]))

expr_sm <- matrix(rep(2,30), ncol=5)
rownames(expr_sm) <- paste('gene', 2:(dim(expr_sm)[1]+1))

expr_sm_projected_back <- projectFromNetworkRecombine(expr, expr_sm)

test_that("projectFromNetworkRecombine uses space of original expression", {
    expect_equal(dim(expr), dim(expr_sm_projected_back))
})

test_that("projectFromNetworkRecombine keeps original values", {
    not_smoothed_genes <- setdiff(rownames(expr), rownames(expr_sm))
    expect_equal(expr_sm_projected_back[not_smoothed_genes,],
                 expr[not_smoothed_genes,])
})

test_that("projectFromNetworkRecombine uses smoothed values", {
    smoothed_genes <- intersect(rownames(expr), rownames(expr_sm))
    expect_equal(expr_sm_projected_back[smoothed_genes,],
                 expr_sm[smoothed_genes,])
})

############################

smoother <- function(expr, adjmatrix, alpha) expr+alpha
expr <- matrix(rep(1,20), ncol=5)
rownames(expr) <- paste('gene', 1:(dim(expr)[1]))

adj <- matrix(rep(1,36), ncol=6)
rownames(adj) <- colnames(adj) <- paste('gene', 2:(dim(adj)[1]+1))

expr_sm <- smoothAndRecombine(expr, adj, .1, smoother)

test_that("smoothAndRecombine smoothes genes in network", {
    smoothed_genes <- intersect(rownames(expr), rownames(adj))
    expect_false(isTRUE(all.equal(expr_sm[smoothed_genes,],
                                  expr[smoothed_genes,])))
})

test_that("smoothAndRecombine leaves unsmoothed genes as is", {
    not_smoothed_genes <- setdiff(rownames(expr), rownames(adj))
    expect_equal(expr_sm[not_smoothed_genes,], expr[not_smoothed_genes,])
})
