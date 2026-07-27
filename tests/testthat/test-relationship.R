be <- simulateBreeding(n_ind = 30, n_marker = 300, seed = 1)

test_that("G is a valid relationship matrix", {
  G <- Gmatrix(be)
  expect_equal(dim(G), c(30L, 30L))
  expect_true(isSymmetric(unname(G), tol = 1e-8))
  expect_equal(rownames(G), colnames(be))
  expect_true(abs(mean(diag(G)) - 1) < 0.4)
})

test_that("G tracks the pedigree relationships it was simulated from", {
  G <- Gmatrix(be)
  A <- Amatrix(be)[colnames(be), colnames(be)]
  expect_gt(cor(A[upper.tri(A)], G[upper.tri(G)]), 0.5)
})

test_that("D is symmetric and differs from G", {
  G <- Gmatrix(be); D <- Dmatrix(be)
  expect_true(isSymmetric(unname(D), tol = 1e-8))
  expect_gt(max(abs(D - G)), 1e-6)
})

test_that("min_maf drops markers and still returns a valid matrix", {
  G <- Gmatrix(be, min_maf = 0.2)
  expect_equal(dim(G), c(30L, 30L))
  expect_true(all(is.finite(G)))
})

test_that("missing genotypes are handled or refused as asked", {
  bem <- simulateBreeding(n_ind = 20, n_marker = 100, missing = 0.05, seed = 2)
  expect_true(all(is.finite(Gmatrix(bem, impute = TRUE))))
  expect_error(Gmatrix(bem, impute = FALSE), "missing")
})

test_that("a monomorphic panel is refused rather than dividing by zero", {
  g <- matrix(1, 10, 5, dimnames = list(NULL, paste0("i", 1:5)))
  expect_error(Gmatrix(BreedingExperiment(g)), "polymorphic")
})
