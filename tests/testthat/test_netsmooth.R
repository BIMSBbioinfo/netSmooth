context("netSmooth")
require(Matrix)

data("smallPPI")
data("smallscRNAseq")

test_that("netSmooth accepts SingleCellExperiment", {
    sink(ifelse(.Platform$OS.type=='unix', '/dev/null', 'NUL'))

    netSmooth(smallscRNAseq, smallPPI, alpha=.5)

    sink()
})

test_that("netSmooth accepts sparse matrices", {
    sink(ifelse(.Platform$OS.type=='unix', '/dev/null', 'NUL'))

    netSmooth(Matrix(assay(smallscRNAseq)), smallPPI, alpha=.5)
    netSmooth(Matrix(assay(smallscRNAseq)), smallPPI, alpha='auto', autoAlphaMethod='entropy')

    sink()
})

test_that("stops with error message if PPI has zero rows/columns", {
    smallPPI.with.zero.row <- data.table::copy(smallPPI)
    smallPPI.with.zero.row[10,] <- 0
    smallPPI.with.zero.row[,10] <- 0
    expect_error(netSmooth(smallscRNAseq, smallPPI.with.zero.row, 0.5))
})


test_that("pickDimReduction is usable with umap", {
  sink(ifelse(.Platform$OS.type=='unix', '/dev/null', 'NUL'))

  expect_equal(pickDimReduction((assay(smallscRNAseq)),
                                flavors='umap',
                                is.counts=TRUE),
               "umap")

  sink()
})

test_that("robusClusters work with umap", {
  sink(ifelse(.Platform$OS.type=='unix', '/dev/null', 'NUL'))

  expect_is(robustClusters(smallscRNAseq,
                           makeConsensusMinSize=2,
                           makeConsensusProportion=.9,
                           dimReduceFlavor = 'umap'),
            "list")

  sink()
})

test_that("netSmooth accepts DelayedMatrix objects", {
  sink(ifelse(.Platform$OS.type=='unix', '/dev/null', 'NUL'))

  singleHDF5 <- as(assay(smallscRNAseq), "HDF5Array")

  rownames(singleHDF5) <- rownames(smallscRNAseq)
  colnames(singleHDF5) <- colnames(smallscRNAseq)

  netSmooth(singleHDF5, smallPPI, alpha=.5)
  netSmooth(singleHDF5, smallPPI, alpha='auto', autoAlphaMethod='entropy')

  sink()
})


