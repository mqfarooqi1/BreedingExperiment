## Assisted-by: Claude (Anthropic). Written with AI assistance under the
## author's direction; methods are established techniques cited in the
## documentation, and the results are validated in tests/testthat.

## Pedigree ordering and the numerator relationship matrix ------------------

#' Order a pedigree so that parents precede their offspring
#'
#' The tabular method for building a relationship matrix requires that every
#' individual appears after both of its parents. This performs that topological
#' sort and reports a cycle if the pedigree contains one.
#'
#' @param ped A data frame or [S4Vectors::DataFrame] with columns `id`, `sire`
#'   and `dam`.
#' @return A data frame with the same columns, reordered.
#' @examples
#' ped <- data.frame(id = c("c", "a", "b"),
#'                   sire = c("a", NA, NA),
#'                   dam = c("b", NA, NA))
#' sortPedigree(ped)
#' @export
sortPedigree <- function(ped) {
    p <- as.data.frame(ped, stringsAsFactors = FALSE)
    for (nm in .PED_COLS) p[[nm]] <- as.character(p[[nm]])
    n <- nrow(p)
    if (n == 0L) return(p)

    placed <- character(0)
    order_idx <- integer(0)
    remaining <- seq_len(n)

    repeat {
        ok <- vapply(remaining, function(i) {
            s <- p$sire[i]; d <- p$dam[i]
            (is.na(s) || s %in% placed) && (is.na(d) || d %in% placed)
        }, logical(1))
        if (!any(ok)) break
        order_idx <- c(order_idx, remaining[ok])
        placed <- c(placed, p$id[remaining[ok]])
        remaining <- remaining[!ok]
        if (!length(remaining)) break
    }
    if (length(remaining)) {
        stop("the pedigree contains a cycle involving: ",
             paste(utils::head(p$id[remaining], 5), collapse = ", "),
             call. = FALSE)
    }
    p[order_idx, , drop = FALSE]
}

#' Pedigree numerator relationship matrix
#'
#' Builds the additive relationship matrix \eqn{A} from the pedigree using the
#' tabular method. Element \eqn{A_{ij}} is twice the coefficient of kinship
#' between individuals \eqn{i} and \eqn{j}, so a diagonal entry is
#' \eqn{1 + F_i} where \eqn{F_i} is the inbreeding coefficient.
#'
#' The matrix is built recursively over a pedigree sorted so that parents come
#' first: an individual's relationship to an earlier individual is the average
#' of its parents' relationships to that individual, and its diagonal is
#' \eqn{1 + \tfrac{1}{2}A_{sire,dam}}.
#'
#' @param x A [BreedingExperiment-class] object.
#' @param ids Optional character vector of identifiers to return, in that
#'   order. Defaults to every individual in the pedigree.
#' @param ... Unused.
#' @return A symmetric numeric matrix with dimnames taken from the pedigree.
#' @references Henderson, C. R. (1976) "A simple method for computing the
#'   inverse of a numerator relationship matrix used in prediction of breeding
#'   values." Biometrics 32, 69-83. \doi{10.2307/2529339}
#' @seealso [Gmatrix()], [Hmatrix()], [inbreeding()]
#' @examples
#' set.seed(1)
#' be <- simulateBreeding(n_ind = 12, n_marker = 30)
#' A <- Amatrix(be)
#' round(A[1:5, 1:5], 3)
#' @name Amatrix
#' @aliases Amatrix Amatrix,BreedingExperiment-method
#' @export
setMethod("Amatrix", "BreedingExperiment", function(x, ids = NULL, ...) {
    ped <- pedigree(x)
    if (nrow(ped) == 0L) {
        stop("this object has no pedigree; see `pedigree<-`.", call. = FALSE)
    }
    A <- .amatrix_from_ped(as.data.frame(ped, stringsAsFactors = FALSE))
    if (is.null(ids)) return(A)
    missing_ids <- setdiff(ids, rownames(A))
    if (length(missing_ids)) {
        stop("not in the pedigree: ",
             paste(utils::head(missing_ids, 5), collapse = ", "), call. = FALSE)
    }
    A[ids, ids, drop = FALSE]
})

#' @keywords internal
#' @noRd
.amatrix_from_ped <- function(p) {
    p <- sortPedigree(p)
    ids <- p$id
    n <- length(ids)
    A <- matrix(0, n, n, dimnames = list(ids, ids))
    s <- match(p$sire, ids)   # NA when unknown
    d <- match(p$dam, ids)

    for (i in seq_len(n)) {
        si <- s[i]; di <- d[i]
        if (i > 1L) {
            j <- seq_len(i - 1L)
            a_s <- if (is.na(si)) 0 else A[j, si]
            a_d <- if (is.na(di)) 0 else A[j, di]
            v <- 0.5 * (a_s + a_d)
            A[j, i] <- v
            A[i, j] <- v
        }
        A[i, i] <- 1 + if (is.na(si) || is.na(di)) 0 else 0.5 * A[si, di]
    }
    A
}

#' Inbreeding coefficients
#'
#' Returns the inbreeding coefficient of each individual, obtained from the
#' diagonal of the pedigree relationship matrix as \eqn{F_i = A_{ii} - 1}.
#'
#' @param x A [BreedingExperiment-class] object.
#' @param ... Unused.
#' @return A named numeric vector of inbreeding coefficients.
#' @seealso [Amatrix()]
#' @examples
#' set.seed(1)
#' be <- simulateBreeding(n_ind = 12, n_marker = 30)
#' round(inbreeding(be), 4)
#' @name inbreeding
#' @aliases inbreeding inbreeding,BreedingExperiment-method
#' @export
setMethod("inbreeding", "BreedingExperiment", function(x, ...) {
    A <- Amatrix(x)
    stats::setNames(diag(A) - 1, rownames(A))
})
