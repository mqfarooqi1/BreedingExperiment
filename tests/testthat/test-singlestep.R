set.seed(1)
be <- simulateBreeding(n_ind = 40, n_marker = 400, genotyped = 25)

test_that("H spans the whole pedigree, not just the genotyped part", {
  H <- Hmatrix(be)
  expect_equal(nrow(H), nrow(pedigree(be)))
  expect_equal(length(attr(H, "genotyped")), 25L)
  expect_true(isSymmetric(unname(H), tol = 1e-6))
})

test_that("H is positive definite", {
  H <- Hmatrix(be)
  expect_gt(min(eigen(H, only.values = TRUE)$values), -1e-6)
})

test_that("blending fully towards the pedigree reproduces A22", {
  H <- Hmatrix(be, blend = 1, tune = FALSE)
  gen <- attr(H, "genotyped")
  A22 <- Amatrix(be)[gen, gen]
  expect_equal(unname(H[gen, gen]), unname(A22), tolerance = 1e-8)
})

test_that("with no blending the genotyped block is the tuned G", {
  H <- Hmatrix(be, blend = 0, tune = FALSE)
  gen <- attr(H, "genotyped")
  expect_equal(unname(H[gen, gen]), unname(Gmatrix(be)[gen, gen]),
               tolerance = 1e-8)
})

test_that("tuning moves G onto the scale of A22", {
  gen <- colnames(be)
  A22 <- Amatrix(be)[gen, gen]
  G <- Gmatrix(be)[gen, gen]
  Gt <- BreedingExperiment:::.tune_G(G, A22)
  off <- function(M) mean(M[upper.tri(M)])
  expect_lt(abs(mean(diag(Gt)) - mean(diag(A22))), 1e-8)
  expect_lt(abs(off(Gt) - off(A22)), 1e-8)
})

test_that("an out-of-range blend is rejected", {
  expect_error(Hmatrix(be, blend = 1.5), "between 0 and 1")
})

test_that("H needs a pedigree", {
  g <- matrix(rbinom(50, 2, .5), 10, 5, dimnames = list(NULL, paste0("i", 1:5)))
  expect_error(Hmatrix(BreedingExperiment(g)), "pedigree")
})

test_that("when everyone is genotyped H reduces to the blended G", {
  full <- simulateBreeding(n_ind = 20, n_marker = 200)
  H <- Hmatrix(full)
  expect_equal(nrow(H), 20L)
})
