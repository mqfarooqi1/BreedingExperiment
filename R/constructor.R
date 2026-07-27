## Assisted-by: Claude (Anthropic). Written with AI assistance under the
## author's direction; methods are established techniques cited in the
## documentation, and the results are validated in tests/testthat.

#' Create a breeding experiment
#'
#' Builds a [BreedingExperiment-class] from a marker genotype matrix, with
#' optional phenotypes, pedigree and marker coordinates.
#'
#' Genotypes are allele dosages: usually the integers 0, 1 and 2, with `NA` for
#' missing calls. Fractional values in between are accepted, so expected dosages
#' from imputation can be stored directly. Markers are rows and individuals are
#' columns, matching the Bioconductor convention that features are rows.
#'
#' If `rowRanges` is not supplied, markers are given placeholder coordinates on
#' a single sequence so that the object is still a valid ranged experiment; pass
#' real coordinates whenever you have them, because that is what allows a data
#' set to be subset by genomic region.
#'
#' The pedigree may name ancestors that were never genotyped. Those extra
#' individuals are exactly what [Hmatrix()] needs in order to combine pedigree
#' and genomic information.
#'
#' @param genotypes A numeric matrix of allele dosages, markers in rows and
#'   individuals in columns, with column names identifying the individuals.
#' @param phenotypes Optional data frame or [S4Vectors::DataFrame] of
#'   individual-level variables, with one row per column of `genotypes`.
#' @param pedigree Optional data frame with columns `id`, `sire` and `dam`.
#'   Unknown parents may be `NA` or `0`.
#' @param rowRanges Optional [GenomicRanges::GRanges-class] giving the position
#'   of each marker, of the same length as `nrow(genotypes)`.
#' @param ... Further assays, passed to
#'   [SummarizedExperiment::SummarizedExperiment()].
#' @return A [BreedingExperiment-class] object.
#' @seealso [Gmatrix()], [Amatrix()], [Hmatrix()], [simulateBreeding()]
#' @examples
#' set.seed(1)
#' geno <- matrix(rbinom(40 * 12, 2, 0.3), nrow = 40, ncol = 12,
#'                dimnames = list(paste0("snp", 1:40), paste0("ind", 1:12)))
#' pheno <- data.frame(yield = rnorm(12))
#' be <- BreedingExperiment(geno, phenotypes = pheno)
#' be
#' dim(be)
#' @export
BreedingExperiment <- function(genotypes, phenotypes = NULL, pedigree = NULL,
                               rowRanges = NULL, ...) {
    genotypes <- .as_geno_matrix(genotypes)
    n_marker <- nrow(genotypes)
    n_ind <- ncol(genotypes)

    if (is.null(colnames(genotypes))) {
        colnames(genotypes) <- paste0("ind", seq_len(n_ind))
    }
    if (is.null(rownames(genotypes))) {
        rownames(genotypes) <- paste0("marker", seq_len(n_marker))
    }

    cd <- if (is.null(phenotypes)) {
        S4Vectors::DataFrame(row.names = colnames(genotypes))
    } else {
        pd <- methods::as(as.data.frame(phenotypes), "DataFrame")
        if (nrow(pd) != n_ind) {
            stop("`phenotypes` must have one row per individual (",
                 n_ind, ").", call. = FALSE)
        }
        rownames(pd) <- colnames(genotypes)
        pd
    }

    rr <- if (is.null(rowRanges)) {
        GenomicRanges::GRanges(
            seqnames = rep("unknown", n_marker),
            ranges = IRanges::IRanges(start = seq_len(n_marker), width = 1L))
    } else {
        if (length(rowRanges) != n_marker) {
            stop("`rowRanges` must have one range per marker (",
                 n_marker, ").", call. = FALSE)
        }
        rowRanges
    }
    names(rr) <- rownames(genotypes)

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = c(list(genotype = genotypes), list(...)),
        rowRanges = rr, colData = cd)

    ped <- .clean_pedigree(pedigree, colnames(genotypes))
    methods::new("BreedingExperiment", se, pedigree = ped)
}

#' @keywords internal
#' @noRd
.as_geno_matrix <- function(g) {
    if (is.data.frame(g)) g <- as.matrix(g)
    if (!is.matrix(g)) stop("`genotypes` must be a matrix.", call. = FALSE)
    storage.mode(g) <- "double"
    ## Dosages are usually the integers 0, 1 and 2, but imputation - whether by
    ## imputeMarkers() or by external software reporting expected dosages -
    ## produces values in between, so anything within [0, 2] is accepted.
    bad <- !is.na(g) & (g < 0 | g > 2)
    if (any(bad)) {
        stop("`genotypes` must be allele dosages between 0 and 2 (or NA).",
             call. = FALSE)
    }
    g
}

#' @keywords internal
#' @noRd
.clean_pedigree <- function(pedigree, ids) {
    if (is.null(pedigree)) {
        return(S4Vectors::DataFrame(id = character(), sire = character(),
                                    dam = character()))
    }
    p <- as.data.frame(pedigree, stringsAsFactors = FALSE)
    if (!all(.PED_COLS %in% names(p))) {
        stop("`pedigree` must have columns 'id', 'sire' and 'dam'.",
             call. = FALSE)
    }
    p <- p[, .PED_COLS, drop = FALSE]
    for (nm in .PED_COLS) p[[nm]] <- as.character(p[[nm]])
    ## treat "0" and "" as unknown, the convention in pedigree files
    p$sire[p$sire %in% c("0", "")] <- NA_character_
    p$dam[p$dam %in% c("0", "")] <- NA_character_

    ## parents that never appear as an id are added as founders, so the
    ## pedigree is self-contained and can be sorted
    parents <- unique(c(p$sire, p$dam))
    parents <- parents[!is.na(parents) & !(parents %in% p$id)]
    if (length(parents)) {
        p <- rbind(data.frame(id = parents, sire = NA_character_,
                              dam = NA_character_, stringsAsFactors = FALSE), p)
    }
    p <- p[!duplicated(p$id), , drop = FALSE]
    S4Vectors::DataFrame(p)
}
