## Assisted-by: Claude (Anthropic). Written with AI assistance under the
## author's direction; methods are established techniques cited in the
## documentation, and the results are validated in tests/testthat.

#' @importFrom methods setClass setGeneric setMethod setValidity new is
#'   validObject callNextMethod slot slotNames show as
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importFrom SummarizedExperiment SummarizedExperiment assay assays
#'   assayNames rowRanges colData rowData
#' @importFrom S4Vectors DataFrame metadata SimpleList
#' @importFrom GenomicRanges GRanges granges
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics ncol nrow
NULL

#' Breeding experiment container
#'
#' A `BreedingExperiment` keeps everything a quantitative-genetic analysis needs
#' in one object: marker genotypes, the phenotypes and design variables of the
#' individuals, the pedigree, and the position of every marker on the genome.
#'
#' The class extends
#' [SummarizedExperiment::RangedSummarizedExperiment-class], so markers are rows
#' and individuals are columns, and the familiar accessors work as usual:
#' `assay()` returns the genotype matrix, `rowRanges()` the marker coordinates
#' as a [GenomicRanges::GRanges-class], and `colData()` the phenotypes. Because
#' marker positions are genomic ranges, a breeding data set can be subset by
#' region, overlapped with annotation, and otherwise handled with the standard
#' Bioconductor vocabulary.
#'
#' One slot is added to the parent class:
#'
#' \describe{
#'   \item{`pedigree`}{A [S4Vectors::DataFrame] with columns `id`, `sire` and
#'     `dam`, giving the parents of each individual. Unknown parents are `NA`.
#'     It may contain ancestors that were never genotyped, which is what makes
#'     single-step analysis possible.}
#' }
#'
#' Genotypes are stored as allele dosages: 0, 1 or 2 copies of the counted
#' allele, with `NA` for missing calls.
#'
#' @section Construction:
#' Use [BreedingExperiment()]; see its help page for examples.
#'
#' @seealso [BreedingExperiment()], [Amatrix()], [Gmatrix()], [Hmatrix()]
#' @name BreedingExperiment-class
#' @aliases BreedingExperiment-class
#' @exportClass BreedingExperiment
setClass("BreedingExperiment",
         contains = "RangedSummarizedExperiment",
         slots = c(pedigree = "DataFrame"),
         prototype = prototype(
           pedigree = S4Vectors::DataFrame(id = character(),
                                           sire = character(),
                                           dam = character())))

.PED_COLS <- c("id", "sire", "dam")

setValidity("BreedingExperiment", function(object) {
  msg <- character()
  ped <- object@pedigree

  if (!all(.PED_COLS %in% colnames(ped))) {
    msg <- c(msg, "`pedigree` must have columns 'id', 'sire' and 'dam'.")
    return(msg)
  }
  if (nrow(ped) > 0L) {
    if (anyDuplicated(ped$id) > 0L) {
      msg <- c(msg, "`pedigree` has duplicated ids.")
    }
    if (anyNA(ped$id)) {
      msg <- c(msg, "`pedigree` has missing ids.")
    }
    known <- c(as.character(ped$sire), as.character(ped$dam))
    known <- known[!is.na(known)]
    if (!all(known %in% as.character(ped$id))) {
      msg <- c(msg,
        "every sire and dam in `pedigree` must itself appear as an id.")
    }
  }

  a <- SummarizedExperiment::assayNames(object)
  if (length(a) == 0L || !"genotype" %in% a) {
    msg <- c(msg, "a 'genotype' assay is required.")
  }
  if (length(msg)) msg else TRUE
})
