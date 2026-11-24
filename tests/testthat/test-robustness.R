## Regression tests for the edge cases raised in Bioconductor review:
## fractional dosages, markers with no calls, silent data handling in the
## constructor, and the simulator's boundaries.

test_that("Dmatrix reproduces the Vitezica codes and interpolates between them", {
    ## one marker, four individuals, including a dosage that is not a whole
    ## number such as mean imputation produces
    g <- matrix(c(0, 1, 2, 1.5), nrow = 1,
                dimnames = list("m1", c("a", "b", "c", "d")))
    be <- BreedingExperiment(g)
    p <- mean(g) / 2
    q <- 1 - p
    code <- function(v) -v^2 + (2 * p + 1) * v - 2 * p^2

    ## exact at the three genotype classes
    expect_equal(code(2), -2 * q^2)
    expect_equal(code(1), 2 * p * q)
    expect_equal(code(0), -2 * p^2)

    ## an expected dosage lands between the classes it sits between, and in
    ## particular is not given the code of the opposite homozygote, which is
    ## what testing dosages for exact equality used to do
    expect_lt(code(1.5), code(1))
    expect_gt(code(1.5), code(2))
    expect_false(isTRUE(all.equal(code(1.5), -2 * p^2)))

    h <- code(as.numeric(g))
    expect_equal(unname(Dmatrix(be)), outer(h, h) / sum((2 * p * q)^2))
})

test_that("Dmatrix gives a missing call no dominance weight", {
    g <- matrix(c(0, 1, 2, 1,
                  2, 1, 0, 1), nrow = 2, byrow = TRUE,
                dimnames = list(c("m1", "m2"), paste0("i", 1:4)))
    gm <- g
    gm[1, 2] <- NA
    be <- BreedingExperiment(gm)
    D <- Dmatrix(be)
    expect_true(all(is.finite(D)))
    expect_true(isSymmetric(unname(D), tol = 1e-8))
    expect_error(Dmatrix(be, impute = FALSE), "missing")
})

test_that("a marker with no calls is dropped rather than invented", {
    set.seed(4)
    g <- matrix(stats::rbinom(10 * 6, 2, 0.4), 10, 6,
                dimnames = list(paste0("m", 1:10), paste0("i", 1:6)))
    g[3, ] <- NA_real_
    be <- BreedingExperiment(g)

    expect_warning(imp <- imputeMarkers(be), "no genotype calls")
    expect_equal(nrow(imp), 9L)
    expect_false("m3" %in% rownames(imp))
    expect_false(anyNA(genotypes(imp)))
    expect_error(imputeMarkers(be, drop_uncalled = FALSE),
                 "nothing to impute")

    ## the relationship matrices must not quietly fill with NaN either
    expect_warning(G <- Gmatrix(be), "no genotype calls")
    expect_true(all(is.finite(G)))
    expect_warning(D <- Dmatrix(be), "no genotype calls")
    expect_true(all(is.finite(D)))
})

test_that("observed heterozygosity counts calls, not imputed dosages", {
    g <- matrix(c(0, 1, 1, 2), nrow = 1,
                dimnames = list("m", paste0("i", 1:4)))
    expect_equal(unname(alleleFrequency(BreedingExperiment(g))$heterozygosity),
                 0.5)

    set.seed(5)
    be <- simulateBreeding(n_ind = 40, n_marker = 30, missing = 0.2)
    imp <- imputeMarkers(be)
    before <- alleleFrequency(be)$heterozygosity
    after <- alleleFrequency(imp)$heterozygosity
    expect_equal(before, after)
    expect_true(all(after >= 0 & after <= 1))

    ## the marker means here land on whole numbers often enough that the
    ## value of an imputed cell cannot be relied on to give it away
    expect_true(any(SummarizedExperiment::assay(imp, "imputed") &
                        genotypes(imp) == round(genotypes(imp))))
})

test_that("imputeMarkers records which cells it filled", {
    set.seed(12)
    be <- simulateBreeding(n_ind = 20, n_marker = 40, missing = 0.15)
    missing_before <- is.na(genotypes(be))
    imp <- imputeMarkers(be)

    expect_true("imputed" %in% SummarizedExperiment::assayNames(imp))
    mask <- SummarizedExperiment::assay(imp, "imputed")
    expect_equal(unname(mask), unname(missing_before))
    expect_identical(dimnames(mask), dimnames(genotypes(imp)))

    ## the mask is part of the container, so subsetting keeps it aligned
    sub <- imp[5:10, 3:8]
    expect_equal(unname(SummarizedExperiment::assay(sub, "imputed")),
                 unname(missing_before[5:10, 3:8]))

    ## imputing an already complete object changes nothing
    expect_identical(genotypes(imputeMarkers(imp)), genotypes(imp))
})

test_that("phenotypes are matched to individuals by row name", {
    set.seed(3)
    g <- matrix(stats::rbinom(20, 2, 0.5), 5, 4,
                dimnames = list(NULL, c("a", "b", "c", "d")))

    ## row names in a different order are a mapping, not a coincidence
    ph <- data.frame(y = c(10, 20, 30, 40),
                     row.names = c("d", "c", "b", "a"))
    be <- BreedingExperiment(g, phenotypes = ph)
    expect_identical(colnames(be), c("a", "b", "c", "d"))
    expect_equal(SummarizedExperiment::colData(be)$y, c(40, 30, 20, 10))

    ## without row names the order given is the order meant
    be2 <- BreedingExperiment(g, phenotypes = data.frame(y = c(1, 2, 3, 4)))
    expect_equal(SummarizedExperiment::colData(be2)$y, c(1, 2, 3, 4))

    ## anything else is refused instead of relabelled
    bad <- data.frame(y = 1:4, row.names = c("a", "b", "c", "zz"))
    expect_error(BreedingExperiment(g, phenotypes = bad), "do not identify")
})

test_that("duplicated individuals are not resolved silently", {
    g <- matrix(1, 5, 2, dimnames = list(NULL, c("x", "y")))

    conflict <- data.frame(id = c("x", "x", "y"), sire = c("y", NA, NA),
                           dam = c(NA, NA, NA))
    expect_error(BreedingExperiment(g, pedigree = conflict),
                 "different parents")

    same <- data.frame(id = c("x", "x", "y"), sire = NA_character_,
                       dam = NA_character_)
    expect_warning(be <- BreedingExperiment(g, pedigree = same),
                   "more than once")
    expect_equal(nrow(pedigree(be)), 2L)

    expect_error(BreedingExperiment(g, pedigree = data.frame(
        id = c("x", NA), sire = NA_character_, dam = NA_character_)),
        "missing id")

    dupcols <- matrix(1, 5, 2, dimnames = list(NULL, c("x", "x")))
    expect_error(BreedingExperiment(dupcols), "duplicated column names")
})

test_that("Hmatrix says so when a genotyped individual is not in the pedigree", {
    set.seed(6)
    be <- simulateBreeding(n_ind = 20, n_marker = 150, n_founder = 4)
    last <- colnames(be)[ncol(be)]
    ped <- as.data.frame(pedigree(be))
    pedigree(be) <- ped[ped$id != last, , drop = FALSE]

    expect_warning(H <- Hmatrix(be), "do not appear in")
    expect_false(last %in% rownames(H))
})

test_that("two founders never give an individual itself as both parents", {
    set.seed(7)
    for (i in seq_len(40)) {
        be <- simulateBreeding(n_ind = 6, n_marker = 10, n_founder = 2)
        ped <- as.data.frame(pedigree(be))
        kids <- ped[!is.na(ped$sire), , drop = FALSE]
        expect_true(all(kids$sire != kids$dam))
        expect_true(all(c(kids$sire, kids$dam) %in% ped$id))
    }
})

test_that("the simulator validates its arguments", {
    expect_error(simulateBreeding(n_ind = 5, h2 = 1.5), "between 0 and 1")
    expect_error(simulateBreeding(n_ind = 5, h2 = -0.1), "between 0 and 1")
    expect_error(simulateBreeding(n_ind = 5, missing = 2), "between 0 and 1")
    expect_error(simulateBreeding(n_ind = 5, missing = NA), "between 0 and 1")
    expect_error(simulateBreeding(n_ind = 0), "positive integer")
    expect_error(simulateBreeding(n_ind = 5, n_marker = -1),
                 "positive integer")
    expect_error(simulateBreeding(n_ind = 5, genotyped = 0),
                 "positive integer")
    expect_error(simulateBreeding(n_ind = 5, genotyped = 9), "cannot exceed")
    expect_error(simulateBreeding(n_ind = 10, n_founder = 1), "two founders")
})

test_that("a population of one individual is allowed", {
    set.seed(8)
    be <- simulateBreeding(n_ind = 1, n_marker = 12)
    expect_equal(dim(be), c(12L, 1L))
    expect_equal(nrow(pedigree(be)), 1L)
    expect_true(all(is.na(pedigree(be)$sire)))
})

test_that("h2 fixes the genetic and residual variances", {
    set.seed(10)
    be <- simulateBreeding(n_ind = 200, n_marker = 100, n_founder = 20,
                           h2 = 0.6)
    cd <- SummarizedExperiment::colData(be)
    expect_equal(stats::var(cd$trueBV), 0.6)
    expect_equal(stats::var(cd$phenotype - cd$trueBV), 0.4, tolerance = 0.25)

    b0 <- simulateBreeding(n_ind = 60, n_marker = 40, n_founder = 10, h2 = 0)
    expect_true(all(SummarizedExperiment::colData(b0)$trueBV == 0))

    b1 <- simulateBreeding(n_ind = 60, n_marker = 40, n_founder = 10, h2 = 1)
    cd1 <- SummarizedExperiment::colData(b1)
    expect_equal(cd1$phenotype, cd1$trueBV)
})

test_that("Dmatrix does not depend on whether imputeMarkers() has been run", {
    ## An imputed cell holds twice the allele frequency, and the dominance
    ## coding returns exactly 2pq there - the heterozygote value - so without
    ## consulting the imputed assay the matrix would shift when a data set was
    ## imputed. Raised in Bioconductor review of 0.99.5.
    set.seed(42)
    be <- simulateBreeding(n_ind = 60, n_marker = 400, missing = 0.05)
    expect_equal(Dmatrix(be), Dmatrix(imputeMarkers(be)))

    ## and the same for the additive matrix, which has always behaved this way
    expect_equal(Gmatrix(be), Gmatrix(imputeMarkers(be)))

    ## the algebraic reason, stated so the test documents itself
    p <- 0.3
    code <- function(g) -g^2 + (2 * p + 1) * g - 2 * p^2
    expect_equal(code(2 * p), 2 * p * (1 - p))

    ## heavier missingness must not drift either
    set.seed(43)
    be2 <- simulateBreeding(n_ind = 40, n_marker = 300, missing = 0.2)
    expect_equal(Dmatrix(be2), Dmatrix(imputeMarkers(be2)))

    ## imputing twice is still the same answer
    expect_equal(Dmatrix(imputeMarkers(be2)),
                 Dmatrix(imputeMarkers(imputeMarkers(be2))))
})