set.seed(1)
ped3 <- data.frame(
  id   = c("A", "B", "C", "D", "E"),
  sire = c(NA,  NA,  "A", NA,  "C"),
  dam  = c(NA,  NA,  "B", NA,  "D"))
g3 <- matrix(1, 10, 5, dimnames = list(paste0("m", 1:10), c("A","B","C","D","E")))
be3 <- BreedingExperiment(g3, pedigree = ped3)

test_that("A matrix reproduces textbook relationships", {
  A <- Amatrix(be3)
  expect_equal(unname(A["A", "B"]), 0)          # unrelated founders
  expect_equal(unname(A["A", "C"]), 0.5)        # parent-offspring
  expect_equal(unname(A["A", "E"]), 0.25)       # grandparent
  expect_true(all(abs(diag(A) - 1) < 1e-12))    # nobody inbred
  expect_true(isSymmetric(unname(A)))
})

test_that("full sibs are related by one half and their offspring is inbred", {
  ped <- data.frame(id = c("S","D","X","Y","Z"),
                    sire = c(NA, NA, "S", "S", "X"),
                    dam  = c(NA, NA, "D", "D", "Y"))
  g <- matrix(1, 8, 5, dimnames = list(paste0("m", 1:8), ped$id))
  A <- Amatrix(BreedingExperiment(g, pedigree = ped))
  expect_equal(unname(A["X", "Y"]), 0.5)
  expect_equal(unname(A["Z", "Z"]), 1.25)       # F = 0.25
})

test_that("inbreeding() equals the diagonal of A minus one", {
  be <- simulateBreeding(n_ind = 15, n_marker = 30)
  expect_equal(unname(inbreeding(be)), unname(diag(Amatrix(be)) - 1))
})

test_that("sortPedigree puts parents before offspring", {
  ped <- data.frame(id = c("c", "a", "b"), sire = c("a", NA, NA),
                    dam = c("b", NA, NA))
  s <- sortPedigree(ped)
  expect_lt(which(s$id == "a"), which(s$id == "c"))
  expect_lt(which(s$id == "b"), which(s$id == "c"))
})

test_that("a cyclic pedigree is rejected", {
  bad <- data.frame(id = c("x", "y"), sire = c("y", "x"), dam = c(NA, NA))
  expect_error(sortPedigree(bad), "cycle")
})

test_that("unknown parents coded 0 become founders", {
  ped <- data.frame(id = c("p", "q"), sire = c("0", "p"), dam = c("0", "0"))
  g <- matrix(1, 5, 2, dimnames = list(NULL, c("p", "q")))
  be <- BreedingExperiment(g, pedigree = ped)
  expect_true(all(is.na(pedigree(be)$sire[pedigree(be)$id == "p"])))
})

test_that("parents absent from the id column are added", {
  ped <- data.frame(id = "kid", sire = "dad", dam = "mum")
  g <- matrix(1, 5, 1, dimnames = list(NULL, "kid"))
  be <- BreedingExperiment(g, pedigree = ped)
  expect_true(all(c("dad", "mum", "kid") %in% pedigree(be)$id))
})
