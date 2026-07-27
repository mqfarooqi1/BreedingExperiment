## Assisted-by: Claude (Anthropic). Written with AI assistance under the
## author's direction; methods are established techniques cited in the
## documentation, and the results are validated in tests/testthat.

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
#' be <- simulateBreeding(n_ind = 25, n_marker = 100)
#' be
#' head(as.data.frame(colData(be)))
#' @export
simulateBreeding <- function(n_ind = 50, n_marker = 200, n_founder = 10,
                             genotyped = NULL, n_qtl = 20, h2 = 0.4,
                             missing = 0) {
    n_founder <- max(2L, min(as.integer(n_founder), n_ind))
    pop <- .simulate_pedigree_genotypes(n_ind, n_marker, n_founder)
    trait <- .simulate_trait(pop$geno, n_qtl, h2)

    keep <- if (is.null(genotyped)) seq_len(n_ind) else
        seq.int(n_ind - min(as.integer(genotyped), n_ind) + 1L, n_ind)

    g <- pop$geno[, keep, drop = FALSE]
    if (missing > 0) {
        n_na <- round(length(g) * missing)
        if (n_na > 0) g[sample(length(g), n_na)] <- NA_real_
    }

    BreedingExperiment(
        genotypes = g,
        phenotypes = data.frame(phenotype = trait$pheno[keep],
                                trueBV = trait$bv[keep]),
        pedigree = pop$pedigree,
        rowRanges = .simple_map(n_marker))
}

## Drop genotypes through a randomly mated pedigree, so markers and pedigree
## tell the same story.
#' @keywords internal
#' @noRd
.simulate_pedigree_genotypes <- function(n_ind, n_marker, n_founder) {
    ids <- sprintf("ind%03d", seq_len(n_ind))
    p <- stats::runif(n_marker, 0.1, 0.9)
    geno <- matrix(0, n_marker, n_ind,
                   dimnames = list(sprintf("snp%04d", seq_len(n_marker)), ids))

    ## founders are unrelated, in Hardy-Weinberg proportions
    for (j in seq_len(n_founder)) {
        geno[, j] <- stats::rbinom(n_marker, 2, p)
    }
    sire <- rep(NA_character_, n_ind)
    dam <- rep(NA_character_, n_ind)

    ## everyone later inherits one allele from each parent
    if (n_ind > n_founder) {
        for (j in seq.int(n_founder + 1L, n_ind)) {
            pool <- seq_len(j - 1L)
            s <- sample(pool, 1L)
            d <- sample(setdiff(pool, s), 1L)
            sire[j] <- ids[s]
            dam[j] <- ids[d]
            geno[, j] <- .gamete(geno[, s]) + .gamete(geno[, d])
        }
    }
    list(geno = geno,
         pedigree = data.frame(id = ids, sire = sire, dam = dam,
                               stringsAsFactors = FALSE))
}

## A trait controlled by a handful of markers, scaled to the target heritability.
#' @keywords internal
#' @noRd
.simulate_trait <- function(geno, n_qtl, h2) {
    n_marker <- nrow(geno)
    qtl <- sort(sample(seq_len(n_marker), min(n_qtl, n_marker)))
    eff <- stats::rnorm(length(qtl))
    bv <- as.numeric(crossprod(geno[qtl, , drop = FALSE], eff))
    bv <- if (stats::sd(bv) > 0) (bv - mean(bv)) / stats::sd(bv) else bv
    ve <- if (h2 > 0) (1 - h2) / h2 else 1
    list(bv = bv, pheno = bv + stats::rnorm(ncol(geno), sd = sqrt(ve)))
}

## Placeholder marker coordinates spread over five sequences.
#' @keywords internal
#' @noRd
.simple_map <- function(n_marker) {
    GenomicRanges::GRanges(
        seqnames = rep(paste0("chr", 1 + (seq_len(n_marker) - 1) %% 5),
                       length.out = n_marker),
        ranges = IRanges::IRanges(start = seq_len(n_marker) * 1000L,
                                  width = 1L))
}

#' @keywords internal
#' @noRd
.gamete <- function(dosage) {
    ## a parent of dosage 0 or 2 always transmits the same allele; a
    ## heterozygote transmits either with probability one half
    ifelse(dosage == 1, stats::rbinom(length(dosage), 1L, 0.5), dosage / 2)
}
