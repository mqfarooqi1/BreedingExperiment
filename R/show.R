#' Display a breeding experiment
#'
#' Prints the size of the object, the assays it holds, the phenotype columns
#' and a summary of the pedigree, including how many individuals are genotyped.
#'
#' @param object A [BreedingExperiment-class] object.
#' @return `object`, invisibly. Called for the side effect of printing.
#' @examples
#' simulateBreeding(n_ind = 10, n_marker = 20, seed = 1)
#' @export
setMethod("show", "BreedingExperiment", function(object) {
    cat("class: BreedingExperiment\n")
    cat("markers:", nrow(object), " individuals genotyped:", ncol(object), "\n")
    cat("assays(", length(SummarizedExperiment::assayNames(object)), "): ",
        paste(SummarizedExperiment::assayNames(object), collapse = " "),
        "\n", sep = "")

    cn <- colnames(SummarizedExperiment::colData(object))
    cat("phenotypes(", length(cn), "): ",
        if (length(cn)) paste(utils::head(cn, 6), collapse = " ") else "none",
        "\n", sep = "")

    sq <- as.character(GenomeInfoDb_seqnames(object))
    cat("sequences(", length(unique(sq)), "): ",
        paste(utils::head(unique(sq), 6), collapse = " "), "\n", sep = "")

    ped <- pedigree(object)
    if (nrow(ped) > 0L) {
        n_gen <- sum(as.character(ped$id) %in% colnames(object))
        n_found <- sum(is.na(ped$sire) & is.na(ped$dam))
        cat("pedigree: ", nrow(ped), " individuals (", n_gen,
            " genotyped, ", nrow(ped) - n_gen, " not; ", n_found,
            " founders)\n", sep = "")
    } else {
        cat("pedigree: none\n")
    }
    invisible(object)
})

#' @keywords internal
#' @noRd
GenomeInfoDb_seqnames <- function(object) {
    rr <- SummarizedExperiment::rowRanges(object)
    if (length(rr) == 0L) return(character())
    as.character(GenomicRanges::seqnames(rr))
}
