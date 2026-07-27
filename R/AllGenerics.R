## Assisted-by: Claude (Anthropic). Written with AI assistance under the
## author's direction; methods are established techniques cited in the
## documentation, and the results are validated in tests/testthat.

#' @rdname pedigree
#' @export
setGeneric("pedigree", function(x, ...) standardGeneric("pedigree"))

#' @rdname pedigree
#' @export
setGeneric("pedigree<-", function(x, value) standardGeneric("pedigree<-"))

#' @rdname genotypes
#' @export
setGeneric("genotypes", function(x, ...) standardGeneric("genotypes"))

#' @rdname genotypes
#' @export
setGeneric("genotypes<-", function(x, value) standardGeneric("genotypes<-"))

#' @rdname alleleFrequency
#' @export
setGeneric("alleleFrequency",
           function(x, ...) standardGeneric("alleleFrequency"))

#' @rdname Amatrix
#' @export
setGeneric("Amatrix", function(x, ...) standardGeneric("Amatrix"))

#' @rdname Gmatrix
#' @export
setGeneric("Gmatrix", function(x, ...) standardGeneric("Gmatrix"))

#' @rdname Dmatrix
#' @export
setGeneric("Dmatrix", function(x, ...) standardGeneric("Dmatrix"))

#' @rdname Hmatrix
#' @export
setGeneric("Hmatrix", function(x, ...) standardGeneric("Hmatrix"))

#' @rdname inbreeding
#' @export
setGeneric("inbreeding", function(x, ...) standardGeneric("inbreeding"))

#' @rdname filterMarkers
#' @export
setGeneric("filterMarkers", function(x, ...) standardGeneric("filterMarkers"))

#' @rdname imputeMarkers
#' @export
setGeneric("imputeMarkers", function(x, ...) standardGeneric("imputeMarkers"))

#' @rdname genotyped
#' @export
setGeneric("genotyped", function(x, ...) standardGeneric("genotyped"))
