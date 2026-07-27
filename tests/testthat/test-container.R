test_that("the constructor validates its inputs", {
  g <- matrix(rbinom(60, 2, .5), 20, 3, dimnames = list(NULL, c("a","b","c")))
  expect_s4_class(BreedingExperiment(g), "BreedingExperiment")
  expect_error(BreedingExperiment(as.list(g)), "matrix")
  expect_error(BreedingExperiment(matrix(5, 4, 2)), "between 0 and 2")
  expect_error(BreedingExperiment(g, phenotypes = data.frame(y = 1:2)),
               "one row per individual")
  expect_error(BreedingExperiment(g, pedigree = data.frame(id = "a")),
               "columns")
})

test_that("fractional dosages from imputation are accepted", {
  g <- matrix(c(0, 0.5, 1.25, 2), 2, 2, dimnames = list(NULL, c("a", "b")))
  expect_s4_class(BreedingExperiment(g), "BreedingExperiment")
})

test_that("it behaves as a RangedSummarizedExperiment", {
  be <- simulateBreeding(n_ind = 12, n_marker = 40, seed = 1)
  expect_s4_class(be, "RangedSummarizedExperiment")
  expect_equal(dim(be), c(40L, 12L))
  expect_s4_class(SummarizedExperiment::rowRanges(be), "GRanges")
  expect_true("phenotype" %in% colnames(SummarizedExperiment::colData(be)))
})

test_that("subsetting keeps the class and the pedigree", {
  be <- simulateBreeding(n_ind = 12, n_marker = 40, seed = 1)
  sub <- be[1:10, 1:5]
  expect_s4_class(sub, "BreedingExperiment")
  expect_equal(dim(sub), c(10L, 5L))
  expect_equal(nrow(pedigree(sub)), nrow(pedigree(be)))
})

test_that("markers can be selected by genomic region", {
  be <- simulateBreeding(n_ind = 10, n_marker = 50, seed = 1)
  chr1 <- be[as.character(GenomicRanges::seqnames(
    SummarizedExperiment::rowRanges(be))) == "chr1", ]
  expect_gt(nrow(chr1), 0)
  expect_s4_class(chr1, "BreedingExperiment")
})

test_that("accessors round-trip", {
  be <- simulateBreeding(n_ind = 10, n_marker = 20, seed = 1)
  g <- genotypes(be)
  genotypes(be) <- g
  expect_equal(genotypes(be), g)
  expect_identical(genotyped(be), colnames(be))

  p <- as.data.frame(pedigree(be))
  pedigree(be) <- p
  expect_equal(nrow(pedigree(be)), nrow(p))
})

test_that("allele frequency summaries are in range", {
  be <- simulateBreeding(n_ind = 25, n_marker = 60, missing = 0.05, seed = 3)
  st <- alleleFrequency(be)
  expect_equal(nrow(st), nrow(be))
  expect_true(all(st$freq >= 0 & st$freq <= 1))
  expect_true(all(st$maf >= 0 & st$maf <= 0.5))
  expect_true(all(st$callRate >= 0 & st$callRate <= 1))
})

test_that("quality control filters and imputes", {
  be <- simulateBreeding(n_ind = 30, n_marker = 120, missing = 0.08, seed = 4)
  expect_lte(nrow(filterMarkers(be, min_maf = 0.05)), nrow(be))
  expect_false(anyNA(genotypes(imputeMarkers(be))))
  expect_error(filterMarkers(be, min_maf = 0.99), "no markers")
})

test_that("show() prints a summary", {
  be <- simulateBreeding(n_ind = 8, n_marker = 15, seed = 1)
  expect_output(show(be), "BreedingExperiment")
  expect_output(show(be), "pedigree")
})

test_that("the simulator honours its arguments", {
  be <- simulateBreeding(n_ind = 30, n_marker = 80, genotyped = 12, seed = 7)
  expect_equal(ncol(be), 12L)
  expect_equal(nrow(pedigree(be)), 30L)
  a <- simulateBreeding(n_ind = 10, n_marker = 20, seed = 9)
  b <- simulateBreeding(n_ind = 10, n_marker = 20, seed = 9)
  expect_equal(genotypes(a), genotypes(b))
})
