context("netSmooth")
require(Matrix)

data("smallPPI")
data("smallscRNAseq")

test_that("netSmooth accepts SingleCellExperiment", {
    sink(ifelse(.Platform$OS.type=='unix', '/dev/null', 'NUL'))

    netSmooth(smallscRNAseq, smallPPI, alpha=.5)
    netSmooth(smallscRNAseq, smallPPI, alpha='auto')

    sink()
})

test_that("netSmooth accepts spase matrices", {
    sink(ifelse(.Platform$OS.type=='unix', '/dev/null', 'NUL'))

    netSmooth(Matrix(assay(smallscRNAseq)), smallPPI, alpha=.5)
    netSmooth(Matrix(assay(smallscRNAseq)), smallPPI, alpha='auto', autoAlphaMethod='entropy')

    sink()
})
