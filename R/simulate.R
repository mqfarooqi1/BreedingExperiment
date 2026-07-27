#' Simulate a small breeding population
#'
#' Generates a pedigreed population with markers and a heritable trait, for
#' examples, tests and teaching. Founders are unrelated; later generations are
#' produced by mating sampled parents, and each offspring inherits one allele
#' from each parent at every marker, so the simulated genotypes and the pedigree
#' agree with each other.
#'
#' Set `genotyped` below `n_ind` to leave some individuals ungenotyped, which is
#' the situation [Hmatrix()] is designed for.
#'
#' The function draws random numbers but does not set the seed itself. Call
#' [set.seed()] beforehand if you need the same population twice.
#'
#' @param n_ind Number of individuals in total.
#' @param n_marker Number of markers.
#' @param n_founder Number of unrelated founders.
#' @param genotyped Number of individuals with genotypes, taken from the most
#'   recent ones. Defaults to all of them.
#' @param n_qtl Number of markers with an effect on the trait.
#' @param h2 Narrow-sense heritability of the simulated trait.
#' @param missing Proportion of genotype calls set to `NA`.
#' @return A [BreedingExperiment-class] object with a `phenotype` column and the
#'   true breeding value in `colData`.
#' @examples
#' set.seed(1)
#' set.seed(1)
#' be <- simulateBreeding(n_ind = 25, n_marker = 100)
#' be
#' head(as.data.frame(colData(be)))
#' @export
simulateBreeding <- function(n_ind = 50, n_marker = 200, n_founder = 10,
                             genotyped = NULL, n_qtl = 20, h2 = 0.4,
                             missing = 0) {
    n_founder <- max(2L, min(as.integer(n_founder), n_ind))
    ids <- sprintf("ind%03d", seq_len(n_ind))

    p <- stats::runif(n_marker, 0.1, 0.9)
    geno <- matrix(0, n_marker, n_ind, dimnames = list(
        sprintf("snp%04d", seq_len(n_marker)), ids))

    ## founders: unrelated, in Hardy-Weinberg proportions
    for (j in seq_len(n_founder)) {
        geno[, j] <- stats::rbinom(n_marker, 2, p)
    }
    sire <- rep(NA_character_, n_ind)
    dam <- rep(NA_character_, n_ind)

    ## later individuals inherit one allele from each parent
    if (n_ind > n_founder) {
        for (j in seq.int(n_founder + 1L, n_ind)) {
            pool <- seq_len(j - 1L)
            s <- sample(pool, 1L)
            d <- sample(setdiff(pool, s), 1L)
            sire[j] <- ids[s]; dam[j] <- ids[d]
            geno[, j] <- .gamete(geno[, s]) + .gamete(geno[, d])
        }
    }

    qtl <- sort(sample(seq_len(n_marker), min(n_qtl, n_marker)))
    eff <- stats::rnorm(length(qtl))
    bv <- as.numeric(crossprod(geno[qtl, , drop = FALSE], eff))
    bv <- if (stats::sd(bv) > 0) (bv - mean(bv)) / stats::sd(bv) else bv
    ve <- if (h2 > 0) (1 - h2) / h2 else 1
    pheno <- bv + stats::rnorm(n_ind, sd = sqrt(ve))

    ped <- data.frame(id = ids, sire = sire, dam = dam,
                      stringsAsFactors = FALSE)

    keep <- if (is.null(genotyped)) seq_len(n_ind) else
        seq.int(n_ind - min(as.integer(genotyped), n_ind) + 1L, n_ind)

    g <- geno[, keep, drop = FALSE]
    if (missing > 0) {
        n_na <- round(length(g) * missing)
        if (n_na > 0) g[sample(length(g), n_na)] <- NA_real_
    }

    rr <- GenomicRanges::GRanges(
        seqnames = rep(paste0("chr", 1 + (seq_len(n_marker) - 1) %% 5),
                       length.out = n_marker),
        ranges = IRanges::IRanges(start = seq_len(n_marker) * 1000L, width = 1L))

    BreedingExperiment(
        genotypes = g,
        phenotypes = data.frame(phenotype = pheno[keep],
                                trueBV = bv[keep]),
        pedigree = ped,
        rowRanges = rr)
}

#' @keywords internal
#' @noRd
.gamete <- function(dosage) {
    ## a parent of dosage 0 or 2 always transmits the same allele; a
    ## heterozygote transmits either with probability one half
    ifelse(dosage == 1, stats::rbinom(length(dosage), 1L, 0.5), dosage / 2)
}

