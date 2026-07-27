test_that("demoBreeding loads and has the documented shape", {
  data(demoBreeding, package = "BreedingExperiment")
  expect_s4_class(demoBreeding, "BreedingExperiment")
  expect_equal(dim(demoBreeding), c(900L, 120L))
  expect_equal(nrow(pedigree(demoBreeding)), 180L)
})

test_that("demoBreeding has ungenotyped ancestors, which is its purpose", {
  data(demoBreeding, package = "BreedingExperiment")
  ung <- setdiff(pedigree(demoBreeding)$id, colnames(demoBreeding))
  expect_equal(length(ung), 60L)
})

test_that("demoBreeding has the documented columns and missingness", {
  data(demoBreeding, package = "BreedingExperiment")
  cd <- SummarizedExperiment::colData(demoBreeding)
  expect_true(all(c("generation", "sex", "yield", "stature",
                    "trueBV_yield", "trueBV_stature") %in% colnames(cd)))
  expect_true(anyNA(cd$yield))          # females only
  expect_false(anyNA(cd$stature))
  expect_true(anyNA(genotypes(demoBreeding)))
})

test_that("the whole workflow runs on demoBreeding", {
  data(demoBreeding, package = "BreedingExperiment")
  be <- imputeMarkers(filterMarkers(demoBreeding, min_maf = 0.05))
  A <- Amatrix(be); G <- Gmatrix(be); H <- Hmatrix(be)
  expect_equal(nrow(A), 180L)
  expect_equal(nrow(G), 120L)
  expect_equal(nrow(H), 180L)
  gen <- colnames(be)
  A22 <- A[gen, gen]
  expect_gt(cor(A22[upper.tri(A22)], G[upper.tri(G)]), 0.5)
  expect_gt(min(eigen(H, only.values = TRUE)$values), -1e-6)
})

test_that("markers span ten chromosomes and can be subset by region", {
  data(demoBreeding, package = "BreedingExperiment")
  sq <- as.character(GenomicRanges::seqnames(
    SummarizedExperiment::rowRanges(demoBreeding)))
  expect_equal(length(unique(sq)), 10L)
  expect_gt(nrow(demoBreeding[sq == "chr1", ]), 0)
})
