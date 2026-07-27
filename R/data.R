#' A demonstration breeding population
#'
#' A small, five-generation breeding population supplied so that every function
#' in the package can be tried on data that behaves like a real programme.
#'
#' The data are **simulated**, not observed. They are built to have the awkward
#' features that real breeding data has and that tidy textbook examples lack:
#'
#' \itemize{
#'   \item **Only part of the population is genotyped.** The two earliest
#'     generations were born before genotyping began, so 60 of the 180
#'     individuals in the pedigree have no markers. This is the situation
#'     [Hmatrix()] exists for.
#'   \item **Unequal family sizes**, because a few sires are used heavily.
#'   \item **A trait recorded on one sex only.** `yield` is missing for males,
#'     as a milk trait would be.
#'   \item **Missing genotype calls**, at about one per cent, as any chip has.
#'   \item **Markers on ten chromosomes** with irregular spacing, so subsetting
#'     by genomic region is meaningful.
#' }
#'
#' Because the data are simulated, the true breeding values are known and are
#' kept in `colData`, which makes the object useful for checking that a
#' prediction method recovers what it should.
#'
#' @format A [BreedingExperiment-class] with 900 markers and 120 genotyped
#'   individuals, and a pedigree of 180. `colData` holds:
#' \describe{
#'   \item{generation}{Generation number, 3 to 5 for genotyped individuals.}
#'   \item{sex}{`"F"` or `"M"`.}
#'   \item{yield}{A moderately heritable trait (\eqn{h^2 \approx 0.35}),
#'     recorded on females only.}
#'   \item{stature}{A more heritable trait (\eqn{h^2 \approx 0.60}), recorded
#'     on everyone.}
#'   \item{trueBV_yield, trueBV_stature}{The simulated true breeding values.}
#' }
#' @source Simulated by `data-raw/make_demo.R` in the package sources.
#' @seealso [simulateBreeding()] to generate your own, [Hmatrix()] for the
#'   analysis this data set is shaped for.
#' @examples
#' data(demoBreeding)
#' demoBreeding
#'
#' # 60 of the 180 individuals in the pedigree were never genotyped
#' nrow(pedigree(demoBreeding)) - ncol(demoBreeding)
#'
#' # the single-step matrix covers all of them
#' H <- Hmatrix(demoBreeding)
#' dim(H)
"demoBreeding"
