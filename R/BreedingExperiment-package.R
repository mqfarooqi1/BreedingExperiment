#' BreedingExperiment: breeding data and relationship matrices
#'
#' Keeps marker genotypes, phenotypes, pedigree and marker coordinates in one
#' object, and computes the relationship matrices that quantitative genetics
#' runs on.
#'
#' The container, [BreedingExperiment-class], extends
#' [SummarizedExperiment::RangedSummarizedExperiment-class], so markers are
#' genomic ranges and the usual Bioconductor accessors apply. The methods that
#' matter are:
#'
#' \describe{
#'   \item{[Amatrix()]}{pedigree numerator relationship matrix}
#'   \item{[Gmatrix()]}{genomic relationship matrix from markers}
#'   \item{[Dmatrix()]}{dominance relationship matrix}
#'   \item{[Hmatrix()]}{single-step matrix combining pedigree and genomic
#'     information, so genotyped and ungenotyped individuals can be analysed
#'     together}
#' }
#'
#' [simulateBreeding()] produces a small pedigreed population for trying these
#' out, and [filterMarkers()] and [imputeMarkers()] cover routine quality
#' control.
#'
#' @section Provenance:
#' A substantial part of this package was written with the assistance of an AI
#' coding assistant (`Assisted-by: Claude, Anthropic`). The author defined the
#' scope, chose the statistical methods, and reviewed, tested and validated the
#' result, and is responsible for the correctness and maintenance of the code.
#'
#' Every method implemented here is an established, published technique in
#' quantitative genetics, cited in the documentation of the function that
#' implements it. The relationship matrices are checked in the test suite
#' against values that follow from theory rather than from a previous run: a
#' parent-offspring relationship of one half, a grandparent relationship of one
#' quarter, full sibs related by one half, the offspring of a full-sib mating
#' inbred at \eqn{F = 0.25}, a positive definite \eqn{H}, and
#' `Hmatrix(blend = 1)` reproducing \eqn{A_{22}} exactly.
#'
#' @keywords internal
"_PACKAGE"
