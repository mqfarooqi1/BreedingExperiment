## Build the `demoBreeding` dataset shipped with the package.
##
## This is SIMULATED data, designed to look like a real breeding programme
## rather than a tidy textbook example:
##   * five generations, with the earliest animals never genotyped, which is
##     what happens when genotyping starts partway through a programme;
##   * unequal family sizes, because a few sires are used heavily;
##   * two traits of different heritability, one of them recorded only on
##     females, so the phenotype column has genuine missingness;
##   * markers spread over ten chromosomes with realistic spacing;
##   * a small amount of missing genotype calls.
##
## Run from the package root:  Rscript data-raw/make_demo.R
suppressMessages(library(BreedingExperiment))
set.seed(20260726)

n_gen      <- 5
n_founder  <- 30
per_gen    <- c(30, 30, 35, 40, 45)   # individuals born in each generation
n_marker   <- 900
n_chr      <- 10
gen_ungeno <- c(1, 2)                  # these generations were never genotyped

ids   <- character(0)
sire  <- character(0)
dam   <- character(0)
gen   <- integer(0)
sex   <- character(0)

## ---- generation 1: unrelated founders --------------------------------------
g1 <- sprintf("A%04d", seq_len(per_gen[1]))
ids <- g1
sire <- rep(NA_character_, length(g1))
dam  <- rep(NA_character_, length(g1))
gen  <- rep(1L, length(g1))
sex  <- rep(c("M", "F"), length.out = length(g1))

## ---- later generations -----------------------------------------------------
counter <- 0L
for (k in 2:n_gen) {
    prev <- ids[gen == (k - 1L)]
    prev_sex <- sex[gen == (k - 1L)]
    sires <- prev[prev_sex == "M"]
    dams  <- prev[prev_sex == "F"]
    if (!length(sires) || !length(dams)) break

    n_k <- per_gen[k]
    ## a few sires are used far more than the rest, as in a real programme
    w <- stats::rgamma(length(sires), shape = 0.7)
    s_k <- sample(sires, n_k, replace = TRUE, prob = w / sum(w))
    d_k <- sample(dams, n_k, replace = TRUE)

    counter <- counter + 1L
    new <- sprintf("%s%04d", LETTERS[k], seq_len(n_k))
    ids  <- c(ids, new)
    sire <- c(sire, s_k)
    dam  <- c(dam, d_k)
    gen  <- c(gen, rep(k, n_k))
    sex  <- c(sex, sample(c("M", "F"), n_k, replace = TRUE))
}

ped <- data.frame(id = ids, sire = sire, dam = dam, stringsAsFactors = FALSE)
n_ind <- nrow(ped)

## ---- drop genotypes through the pedigree -----------------------------------
p0 <- stats::runif(n_marker, 0.05, 0.95)
G <- matrix(0L, n_marker, n_ind, dimnames = list(NULL, ids))
gamete <- function(d) ifelse(d == 1L, stats::rbinom(length(d), 1L, 0.5), d %/% 2L)

for (i in seq_len(n_ind)) {
    s <- ped$sire[i]; d <- ped$dam[i]
    if (is.na(s) || is.na(d)) {
        G[, i] <- stats::rbinom(n_marker, 2, p0)
    } else {
        G[, i] <- gamete(G[, s]) + gamete(G[, d])
    }
}

## ---- two traits ------------------------------------------------------------
qtl1 <- sort(sample(n_marker, 60))
qtl2 <- sort(sample(n_marker, 25))
eff1 <- stats::rnorm(length(qtl1))
eff2 <- stats::rnorm(length(qtl2), sd = 1.4)

scale_bv <- function(v) if (stats::sd(v) > 0) (v - mean(v)) / stats::sd(v) else v
bv1 <- scale_bv(as.numeric(crossprod(G[qtl1, , drop = FALSE], eff1)))
bv2 <- scale_bv(as.numeric(crossprod(G[qtl2, , drop = FALSE], eff2)))

h2_1 <- 0.35; h2_2 <- 0.60
yield <- 6000 + 400 * bv1 + stats::rnorm(n_ind, sd = 400 * sqrt((1 - h2_1) / h2_1))
stature <- 145 + 5 * bv2 + stats::rnorm(n_ind, sd = 5 * sqrt((1 - h2_2) / h2_2))

## yield is recorded on females only, as for a milk trait
yield[sex == "M"] <- NA_real_

pheno <- data.frame(
    generation = gen,
    sex        = factor(sex, levels = c("F", "M")),
    yield      = round(yield, 1),
    stature    = round(stature, 1),
    trueBV_yield   = round(bv1, 4),
    trueBV_stature = round(bv2, 4),
    stringsAsFactors = FALSE)

## ---- who was genotyped -----------------------------------------------------
keep <- which(!(gen %in% gen_ungeno))
Gk <- G[, keep, drop = FALSE]
storage.mode(Gk) <- "integer"

## a little missingness, as any real chip has
n_na <- round(length(Gk) * 0.01)
Gk[sample(length(Gk), n_na)] <- NA_integer_

## ---- marker map ------------------------------------------------------------
chr <- rep(seq_len(n_chr), length.out = n_marker)
chr <- sort(chr)
pos <- unlist(lapply(split(seq_len(n_marker), chr), function(idx) {
    cumsum(sample(20000:900000, length(idx), replace = TRUE))
}), use.names = FALSE)
rr <- GenomicRanges::GRanges(
    seqnames = paste0("chr", chr),
    ranges = IRanges::IRanges(start = pos, width = 1L))
names(rr) <- sprintf("snp%05d", seq_len(n_marker))
rownames(Gk) <- names(rr)

demoBreeding <- BreedingExperiment(
    genotypes  = Gk,
    phenotypes = pheno[keep, , drop = FALSE],
    pedigree   = ped,
    rowRanges  = rr)

dir.create("data", showWarnings = FALSE)
save(demoBreeding, file = "data/demoBreeding.rda", compress = "xz")

cat(sprintf("markers %d | genotyped %d | pedigree %d | size %.0f KB\n",
            nrow(demoBreeding), ncol(demoBreeding),
            nrow(pedigree(demoBreeding)),
            file.size("data/demoBreeding.rda") / 1024))
