## Genomic relationship matrices ---------------------------------------------

#' Genomic relationship matrix
#'
#' Builds the additive genomic relationship matrix \eqn{G} from marker dosages
#' using VanRaden's first method: with \eqn{W} the marker matrix centred by
#' twice the allele frequency, \eqn{G = WW' / 2\sum p(1-p)}.
#'
#' \eqn{G} plays the same role as the pedigree matrix \eqn{A} but is estimated
#' from markers, so it captures Mendelian sampling that a pedigree cannot: full
#' sibs have an expected pedigree relationship of exactly 0.5, whereas their
#' genomic relationship varies around it.
#'
#' @param x A [BreedingExperiment-class] object.
#' @param min_maf Markers with a minor allele frequency below this are dropped.
#' @param impute Replace missing dosages by twice the allele frequency before
#'   computing, which gives a missing call no weight in the result. If `FALSE`,
#'   missing values raise an error.
#' @param ... Unused.
#' @return A symmetric numeric matrix over the genotyped individuals.
#' @references VanRaden, P. M. (2008) "Efficient methods to compute genomic
#'   predictions." Journal of Dairy Science 91, 4414-4423.
#'   \doi{10.3168/jds.2007-0980}
#' @seealso [Amatrix()], [Dmatrix()], [Hmatrix()]
#' @examples
#' set.seed(1)
#' be <- simulateBreeding(n_ind = 20, n_marker = 200)
#' G <- Gmatrix(be)
#' round(G[1:4, 1:4], 3)
#' @name Gmatrix
#' @aliases Gmatrix Gmatrix,BreedingExperiment-method
#' @export
setMethod("Gmatrix", "BreedingExperiment",
          function(x, min_maf = 0, impute = TRUE, ...) {
    W <- .centred_markers(x, min_maf = min_maf, impute = impute)
    denom <- attr(W, "denom")
    G <- crossprod(W) / denom          # individuals are columns of the assay
    dimnames(G) <- list(colnames(x), colnames(x))
    G
})

#' Dominance genomic relationship matrix
#'
#' Builds the dominance relationship matrix from marker dosages using the
#' genotypic parameterisation of Vitezica, Varona and Legarra (2013).
#' Heterozygotes and the two homozygotes are coded so that the resulting matrix
#' captures dominance deviations, orthogonal to the additive relationships
#' returned by [Gmatrix()].
#'
#' Fitting additive and dominance effects together lets a model separate the
#' part of the genetic variance that is inherited predictably from the part
#' that is not, which matters for crossbreeding and for predicting the
#' performance of specific crosses.
#'
#' @section Dosages that are not whole numbers:
#' Vitezica and others code the genotypes AA, Aa and aa as \eqn{-2q^2},
#' \eqn{2pq} and \eqn{-2p^2}. Those three values lie on a single quadratic in
#' the dosage \eqn{g},
#' \deqn{-g^2 + (2p + 1)g - 2p^2,}
#' which is what this function evaluates. At dosages 2, 1 and 0 it returns the
#' three published values exactly, and an expected dosage from imputation
#' falls between them rather than being forced into one genotype class.
#'
#' A missing call is given a dominance contribution of zero, which is the
#' weight that mean imputation already gives it in [Gmatrix()]. A cell that
#' [imputeMarkers()] has filled is treated the same way, because an expected
#' dosage is not a genotype: at a dosage of \eqn{2p} the quadratic returns
#' exactly \eqn{2pq}, so counting it would record every imputed cell as a
#' heterozygote. The result is therefore the same whether or not
#' [imputeMarkers()] has been run, provided the object carries the `imputed`
#' assay that records which cells were filled.
#'
#' @inheritParams Gmatrix
#' @param impute Allow missing dosages, giving them a dominance contribution
#'   of zero. If `FALSE`, missing values raise an error.
#' @return A symmetric numeric matrix over the genotyped individuals.
#' @references Vitezica, Z. G., Varona, L. & Legarra, A. (2013) "On the
#'   additive and dominant variance and covariance of individuals within the
#'   genomic selection scope." Genetics 195, 1223-1230.
#'   \doi{10.1534/genetics.113.155176}
#' @seealso [Gmatrix()]
#' @examples
#' set.seed(1)
#' be <- simulateBreeding(n_ind = 20, n_marker = 200)
#' D <- Dmatrix(be)
#' round(D[1:4, 1:4], 3)
#' @name Dmatrix
#' @aliases Dmatrix Dmatrix,BreedingExperiment-method
#' @export
setMethod("Dmatrix", "BreedingExperiment",
          function(x, min_maf = 0, impute = TRUE, ...) {
    ## missing calls are kept as NA here rather than filled with 2p, because
    ## the dominance code of an unknown genotype is its expectation, zero,
    ## not the code of the average dosage
    g <- .prepared_genotypes(x, min_maf = min_maf, impute = impute,
                             keep_na = TRUE)
    p <- rowMeans(g, na.rm = TRUE) / 2
    q <- 1 - p
    ## the quadratic through the three Vitezica genotype codes; see Details
    H <- -g^2 + (2 * p + 1) * g - 2 * p^2
    ## A cell with no genotype call carries no information about dominance,
    ## and that is equally true whether the call is still missing or has been
    ## replaced by an expected dosage. Both get a contribution of zero, so
    ## running imputeMarkers() first cannot change the result. Without the
    ## second step an imputed dosage of 2p would be read as a heterozygote,
    ## because the coding above returns exactly 2pq there.
    H[is.na(H)] <- 0
    filled <- .filled_cells(x, g)
    if (!is.null(filled)) H[filled] <- 0
    denom <- sum((2 * p * q)^2)
    if (denom <= 0) {
        stop("all markers are monomorphic; cannot build D.", call. = FALSE)
    }
    D <- crossprod(H) / denom
    dimnames(D) <- list(colnames(x), colnames(x))
    D
})

## Shared genotype preparation: drop markers nothing can be learned from,
## deal with missing calls, then filter on frequency and variation.
##
## `keep_na = TRUE` leaves missing calls in place for callers that give them
## their own weight; every frequency and variance below is computed with
## `na.rm = TRUE` so that the filtering is the same either way.
#' @keywords internal
#' @noRd
.prepared_genotypes <- function(x, min_maf = 0, impute = TRUE,
                                keep_na = FALSE) {
    g <- .drop_uncalled(genotypes(x))
    if (anyNA(g)) {
        if (!impute) {
            stop("genotypes contain missing values; use impute = TRUE or ",
                 "see `imputeMarkers()`.", call. = FALSE)
        }
        if (!keep_na) {
            mu <- rowMeans(g, na.rm = TRUE)
            idx <- which(is.na(g), arr.ind = TRUE)
            g[idx] <- mu[idx[, "row"]]
        }
    }
    keep <- rep(TRUE, nrow(g))
    if (min_maf > 0) {
        p <- rowMeans(g, na.rm = TRUE) / 2
        keep <- keep & pmin(p, 1 - p) >= min_maf
    }
    ## a marker with no variation contributes nothing and breaks the scaling
    sds <- apply(g, 1L, stats::sd, na.rm = TRUE)
    keep <- keep & !is.na(sds) & sds > 0
    g <- g[keep, , drop = FALSE]
    if (nrow(g) == 0L) {
        stop("no polymorphic markers left after filtering.", call. = FALSE)
    }
    g
}

## Cells that hold an expected dosage rather than a call, aligned to the
## markers that survived preparation. Both the uncalled-marker drop and the
## frequency filters remove rows, so the mask is matched up by marker name.
#' @keywords internal
#' @noRd
.filled_cells <- function(x, g) {
    m <- .imputed_mask(x)
    if (identical(m, FALSE)) return(NULL)
    rn <- rownames(g)
    if (is.null(rn) || is.null(rownames(m))) return(NULL)
    m[rn, , drop = FALSE]
}

## A marker with no calls at all has no allele frequency to impute from, and
## carrying it further turns the relationship matrix into NaN without saying
## why.
#' @keywords internal
#' @noRd
.drop_uncalled <- function(g) {
    uncalled <- rowSums(!is.na(g)) == 0L
    if (any(uncalled)) {
        warning(sum(uncalled), " marker(s) have no genotype calls at all ",
                "and were dropped: ", .some(rownames(g)[uncalled]),
                call. = FALSE)
        g <- g[!uncalled, , drop = FALSE]
    }
    if (nrow(g) == 0L) {
        stop("no marker has any genotype call.", call. = FALSE)
    }
    g
}

#' @keywords internal
#' @noRd
.centred_markers <- function(x, min_maf = 0, impute = TRUE) {
    g <- .prepared_genotypes(x, min_maf = min_maf, impute = impute)
    p <- rowMeans(g) / 2
    W <- g - 2 * p                       # recycles over rows (markers)
    denom <- 2 * sum(p * (1 - p))
    if (denom <= 0) {
        stop("all markers are monomorphic; cannot build G.", call. = FALSE)
    }
    attr(W, "denom") <- denom
    W
}