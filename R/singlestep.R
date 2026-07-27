## Single-step relationships: combining pedigree and genomic information -----

#' Single-step relationship matrix
#'
#' Builds the matrix \eqn{H} that combines pedigree and genomic information, so
#' that genotyped and ungenotyped individuals can be analysed in a single
#' model. This is the relationship structure behind single-step genomic
#' prediction (ssGBLUP).
#'
#' Real breeding programmes genotype only part of the population. Using only
#' \eqn{G} throws away every ungenotyped relative; using only \eqn{A} throws
#' away the markers. \eqn{H} keeps both: relationships among genotyped
#' individuals come from the markers, relationships among ungenotyped
#' individuals are those from the pedigree corrected by what the markers reveal
#' about their genotyped relatives, and the two groups are connected through the
#' pedigree.
#'
#' Writing 1 for ungenotyped and 2 for genotyped individuals,
#' \deqn{H_{22} = G^{*}}
#' \deqn{H_{12} = A_{12} A_{22}^{-1} G^{*}}
#' \deqn{H_{11} = A_{11} + A_{12} A_{22}^{-1} (G^{*} - A_{22}) A_{22}^{-1} A_{21}}
#'
#' Because \eqn{G} and \eqn{A_{22}} are estimated on different scales, \eqn{G}
#' is first made compatible with \eqn{A_{22}} in two standard steps. *Tuning*
#' rescales \eqn{G} so that its mean diagonal and mean off-diagonal match those
#' of \eqn{A_{22}}. *Blending* then mixes in a small share of \eqn{A_{22}},
#' \eqn{G^{*} = (1-w)G + wA_{22}}, which also makes the result invertible when
#' markers are fewer than individuals.
#'
#' @param x A [BreedingExperiment-class] object with both a pedigree and
#'   genotypes.
#' @param blend Weight \eqn{w} given to \eqn{A_{22}} when blending, between 0
#'   and 1. The default, 0.05, is a common choice.
#' @param tune Rescale \eqn{G} to the scale of \eqn{A_{22}} before blending.
#' @param min_maf,impute Passed to [Gmatrix()].
#' @param ... Unused.
#' @return A symmetric numeric matrix over every individual in the pedigree,
#'   ungenotyped individuals first, with an attribute `genotyped` giving the
#'   identifiers of the genotyped group.
#' @references
#' Legarra, A., Aguilar, I. & Misztal, I. (2009) "A relationship matrix
#' including full pedigree and genomic information." Journal of Dairy Science
#' 92, 4656-4663. \doi{10.3168/jds.2009-2061}
#'
#' Christensen, O. F. & Lund, M. S. (2010) "Genomic prediction when some
#' animals are not genotyped." Genetics Selection Evolution 42, 2.
#' \doi{10.1186/1297-9686-42-2}
#'
#' Aguilar, I., Misztal, I., Johnson, D. L., Legarra, A., Tsuruta, S. &
#' Lawlor, T. J. (2010) "Hot topic: A unified approach to utilize phenotypic,
#' full pedigree, and genomic information for genetic evaluation of Holstein
#' final score." Journal of Dairy Science 93, 743-752.
#' \doi{10.3168/jds.2009-2730}
#' @seealso [Amatrix()], [Gmatrix()]
#' @examples
#' # only some of the individuals are genotyped
#' set.seed(1)
#' be <- simulateBreeding(n_ind = 30, n_marker = 300, genotyped = 20)
#' H <- Hmatrix(be)
#' dim(H)
#' attr(H, "genotyped")
#' @name Hmatrix
#' @aliases Hmatrix Hmatrix,BreedingExperiment-method
#' @export
setMethod("Hmatrix", "BreedingExperiment",
          function(x, blend = 0.05, tune = TRUE, min_maf = 0, impute = TRUE,
                   ...) {
    if (blend < 0 || blend > 1) {
        stop("`blend` must be between 0 and 1.", call. = FALSE)
    }
    ped <- pedigree(x)
    if (nrow(ped) == 0L) {
        stop("a pedigree is required for a single-step matrix.", call. = FALSE)
    }
    A <- Amatrix(x)
    gen <- intersect(rownames(A), colnames(x))
    if (!length(gen)) {
        stop("no genotyped individual appears in the pedigree.", call. = FALSE)
    }
    ung <- setdiff(rownames(A), gen)

    G <- Gmatrix(x, min_maf = min_maf, impute = impute)[gen, gen, drop = FALSE]
    A22 <- A[gen, gen, drop = FALSE]

    if (tune) G <- .tune_G(G, A22)
    Gs <- (1 - blend) * G + blend * A22

    if (!length(ung)) {
        out <- Gs
        attr(out, "genotyped") <- gen
        return(out)
    }

    A11 <- A[ung, ung, drop = FALSE]
    A12 <- A[ung, gen, drop = FALSE]
    A22inv <- .safe_solve(A22, "A22")
    P <- A12 %*% A22inv                      # A12 A22^-1

    H11 <- A11 + P %*% (Gs - A22) %*% t(P)
    H12 <- P %*% Gs
    H <- rbind(cbind(H11, H12), cbind(t(H12), Gs))
    dimnames(H) <- list(c(ung, gen), c(ung, gen))
    H <- (H + t(H)) / 2                      # remove any asymmetry from rounding
    attr(H, "genotyped") <- gen
    H
})

## Rescale G so its mean diagonal and mean off-diagonal match those of A22.
## Following the common two-parameter tuning, G* = a + b G.
#' @keywords internal
#' @noRd
.tune_G <- function(G, A22) {
    n <- nrow(G)
    if (n < 2L) return(G)
    off <- function(M) mean(M[upper.tri(M)])
    dg <- function(M) mean(diag(M))
    mG_d <- dg(G); mG_o <- off(G)
    mA_d <- dg(A22); mA_o <- off(A22)
    denom <- mG_d - mG_o
    if (!is.finite(denom) || abs(denom) < 1e-12) return(G)
    b <- (mA_d - mA_o) / denom
    a <- mA_o - b * mG_o
    b * G + a
}

#' @keywords internal
#' @noRd
.safe_solve <- function(M, what = "matrix") {
    out <- try(solve(M), silent = TRUE)
    if (!inherits(out, "try-error")) return(out)
    ## a pedigree with repeated identical rows can be singular; nudge the
    ## diagonal, which is standard practice, and say so
    eps <- 1e-6 * mean(diag(M))
    out <- try(solve(M + diag(eps, nrow(M))), silent = TRUE)
    if (inherits(out, "try-error")) {
        stop(what, " is singular and could not be inverted.", call. = FALSE)
    }
    warning(what, " was singular; a small ridge was added to invert it.",
            call. = FALSE)
    out
}
