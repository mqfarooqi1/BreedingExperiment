#' Pedigree of a breeding experiment
#'
#' Gets or sets the pedigree, a table with columns `id`, `sire` and `dam`. The
#' pedigree may include ancestors that were never genotyped; those individuals
#' are what allow [Hmatrix()] to combine pedigree and genomic information.
#'
#' @param x A [BreedingExperiment-class] object.
#' @param value A data frame with columns `id`, `sire` and `dam`.
#' @param ... Unused.
#' @return `pedigree()` returns a [S4Vectors::DataFrame]; the replacement form
#'   returns the updated object.
#' @examples
#' be <- simulateBreeding(n_ind = 8, n_marker = 20, seed = 1)
#' head(pedigree(be))
#' @name pedigree
#' @aliases pedigree pedigree<-
NULL

#' @rdname pedigree
#' @export
setMethod("pedigree", "BreedingExperiment", function(x, ...) x@pedigree)

#' @rdname pedigree
#' @export
setReplaceMethod("pedigree", "BreedingExperiment", function(x, value) {
    x@pedigree <- .clean_pedigree(value, colnames(x))
    methods::validObject(x)
    x
})

#' Genotype dosages
#'
#' Gets or sets the matrix of allele dosages, a shorthand for the `genotype`
#' assay. Markers are rows and individuals are columns.
#'
#' @param x A [BreedingExperiment-class] object.
#' @param value A numeric matrix of dosages coded 0, 1 or 2.
#' @param ... Unused.
#' @return `genotypes()` returns a numeric matrix; the replacement form returns
#'   the updated object.
#' @examples
#' be <- simulateBreeding(n_ind = 6, n_marker = 10, seed = 1)
#' genotypes(be)[1:5, 1:4]
#' @name genotypes
#' @aliases genotypes genotypes<-
NULL

#' @rdname genotypes
#' @export
setMethod("genotypes", "BreedingExperiment", function(x, ...) {
    SummarizedExperiment::assay(x, "genotype")
})

#' @rdname genotypes
#' @export
setReplaceMethod("genotypes", "BreedingExperiment", function(x, value) {
    value <- .as_geno_matrix(value)
    SummarizedExperiment::assay(x, "genotype") <- value
    x
})

#' Which individuals are genotyped
#'
#' Returns the identifiers that appear in the genotype matrix. In a single-step
#' analysis the pedigree usually names more individuals than this, and the
#' difference is what [Hmatrix()] reconciles.
#'
#' @param x A [BreedingExperiment-class] object.
#' @param ... Unused.
#' @return A character vector of individual identifiers.
#' @examples
#' be <- simulateBreeding(n_ind = 8, n_marker = 20, seed = 1)
#' genotyped(be)
#' @name genotyped
#' @aliases genotyped genotyped,BreedingExperiment-method
#' @export
setMethod("genotyped", "BreedingExperiment", function(x, ...) colnames(x))

#' Allele frequencies and marker summaries
#'
#' Computes, for every marker, the frequency of the counted allele, the minor
#' allele frequency, the call rate and the observed heterozygosity.
#'
#' @param x A [BreedingExperiment-class] object.
#' @param ... Unused.
#' @return A [S4Vectors::DataFrame] with one row per marker and columns `freq`,
#'   `maf`, `callRate` and `heterozygosity`.
#' @examples
#' be <- simulateBreeding(n_ind = 30, n_marker = 50, seed = 1)
#' head(alleleFrequency(be))
#' @name alleleFrequency
#' @aliases alleleFrequency alleleFrequency,BreedingExperiment-method
#' @export
setMethod("alleleFrequency", "BreedingExperiment", function(x, ...) {
    g <- genotypes(x)
    p <- rowMeans(g, na.rm = TRUE) / 2
    S4Vectors::DataFrame(
        freq = p,
        maf = pmin(p, 1 - p),
        callRate = rowMeans(!is.na(g)),
        heterozygosity = rowMeans(g == 1, na.rm = TRUE),
        row.names = rownames(g))
})
