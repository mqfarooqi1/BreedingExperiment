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
#' @keywords internal
"_PACKAGE"
