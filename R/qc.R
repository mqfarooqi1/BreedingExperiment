#' Filter markers on quality
#'
#' Removes markers that are uninformative or poorly genotyped: those below a
#' minor allele frequency threshold, those called in too few individuals, and
#' optionally those that are monomorphic.
#'
#' @param x A [BreedingExperiment-class] object.
#' @param min_maf Minimum minor allele frequency to keep.
#' @param min_call_rate Minimum proportion of individuals with a call.
#' @param drop_monomorphic Drop markers with no variation.
#' @param ... Unused.
#' @return A [BreedingExperiment-class] object with fewer rows.
#' @seealso [imputeMarkers()], [alleleFrequency()]
#' @examples
#' set.seed(1)
#' be <- simulateBreeding(n_ind = 30, n_marker = 100, missing = 0.05)
#' nrow(be)
#' nrow(filterMarkers(be, min_maf = 0.05, min_call_rate = 0.9))
#' @name filterMarkers
#' @aliases filterMarkers filterMarkers,BreedingExperiment-method
#' @export
setMethod("filterMarkers", "BreedingExperiment",
          function(x, min_maf = 0.01, min_call_rate = 0.9,
                   drop_monomorphic = TRUE, ...) {
    st <- alleleFrequency(x)
    keep <- st$maf >= min_maf & st$callRate >= min_call_rate
    if (drop_monomorphic) {
        g <- genotypes(x)
        varies <- apply(g, 1L, function(v) {
            v <- v[!is.na(v)]
            length(v) > 0L && stats::sd(v) > 0
        })
        keep <- keep & varies
    }
    keep[is.na(keep)] <- FALSE
    if (!any(keep)) stop("no markers passed the filter.", call. = FALSE)
    x[keep, ]
})

#' Impute missing genotypes
#'
#' Fills missing dosages with twice the allele frequency of the marker, the
#' simple mean imputation that relationship matrices assume. It is adequate for
#' the low rates of missingness left after quality control; for anything more,
#' impute with dedicated software before building the object.
#'
#' @param x A [BreedingExperiment-class] object.
#' @param ... Unused.
#' @return A [BreedingExperiment-class] object with no missing genotypes.
#' @seealso [filterMarkers()]
#' @examples
#' set.seed(1)
#' be <- simulateBreeding(n_ind = 20, n_marker = 50, missing = 0.1)
#' anyNA(genotypes(be))
#' anyNA(genotypes(imputeMarkers(be)))
#' @name imputeMarkers
#' @aliases imputeMarkers imputeMarkers,BreedingExperiment-method
#' @export
setMethod("imputeMarkers", "BreedingExperiment", function(x, ...) {
    g <- genotypes(x)
    if (!anyNA(g)) return(x)
    mu <- rowMeans(g, na.rm = TRUE)
    mu[is.nan(mu)] <- 0
    idx <- which(is.na(g), arr.ind = TRUE)
    g[idx] <- mu[idx[, "row"]]
    genotypes(x) <- g
    x
})
