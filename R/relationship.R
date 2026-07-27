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
#'   computing. If `FALSE`, missing values raise an error.
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
#' @inheritParams Gmatrix
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
    g <- .prepared_genotypes(x, min_maf = min_maf, impute = impute)
    p <- rowMeans(g) / 2
    q <- 1 - p
    ## AA -> -2q^2 ; Aa -> 2pq ; aa -> -2p^2   (markers in rows)
    H <- matrix(0, nrow(g), ncol(g), dimnames = dimnames(g))
    H[] <- ifelse(g == 2, -2 * q^2, ifelse(g == 1, 2 * p * q, -2 * p^2))
    denom <- sum((2 * p * q)^2)
    if (denom <= 0) {
        stop("all markers are monomorphic; cannot build D.", call. = FALSE)
    }
    D <- crossprod(H) / denom
    dimnames(D) <- list(colnames(x), colnames(x))
    D
})

#' @keywords internal
#' @noRd
.prepared_genotypes <- function(x, min_maf = 0, impute = TRUE) {
    g <- genotypes(x)
    if (anyNA(g)) {
        if (!impute) {
            stop("genotypes contain missing values; use impute = TRUE or ",
                 "see `imputeMarkers()`.", call. = FALSE)
        }
        mu <- rowMeans(g, na.rm = TRUE)
        idx <- which(is.na(g), arr.ind = TRUE)
        if (nrow(idx)) g[idx] <- mu[idx[, "row"]]
    }
    if (min_maf > 0) {
        p <- rowMeans(g) / 2
        keep <- pmin(p, 1 - p) >= min_maf
        g <- g[keep, , drop = FALSE]
    }
    ## a marker with no variation contributes nothing and breaks the scaling
    keep <- apply(g, 1L, function(v) stats::sd(v) > 0)
    g <- g[keep, , drop = FALSE]
    if (nrow(g) == 0L) {
        stop("no polymorphic markers left after filtering.", call. = FALSE)
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
