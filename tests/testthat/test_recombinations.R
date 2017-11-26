context("Recombinations")

# projectOnNetwork(gene_expression, new_features, missing.value=0)

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
