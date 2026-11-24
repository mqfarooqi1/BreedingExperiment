# Changelog

## BreedingExperiment 0.99.6

### Bug fixes

- [`Dmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Dmatrix.md)
  no longer changes when
  [`imputeMarkers()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/imputeMarkers.md)
  has been run first. The dominance coding returns exactly `2pq` at a
  dosage of `2p`, which is the heterozygote value, so a mean-imputed
  cell was being counted as an observed heterozygote while a
  still-missing call correctly contributed zero. On a data set with 5%
  missing genotypes this shifted entries by up to 0.08 and raised the
  mean diagonal from 0.91 to 0.96. Cells recorded in the `imputed` assay
  now get the same zero contribution as a missing call, so the matrix is
  identical either way. Raised in the Bioconductor review of 0.99.5;
  [`Gmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Gmatrix.md)
  was never affected, because mean imputation already centres an imputed
  cell to zero there. \# BreedingExperiment 0.99.5

Changes made in response to the Bioconductor review
([\#117](https://github.com/Bioconductor/BiocContributions/issues/117)).

### Bug fixes

- [`Dmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Dmatrix.md)
  no longer mis-codes dosages that are not whole numbers. The Vitezica
  genotype codes are now evaluated as the quadratic that passes through
  all three of them, so dosages of 0, 1 and 2 give exactly the published
  values and an expected dosage from imputation falls between them.
  Previously anything other than exactly 1 or 2 was coded as the
  opposite homozygote, which made the dominance matrix wrong whenever
  genotypes were imputed - including under the default `impute = TRUE`.
- A missing call now contributes zero to
  [`Dmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Dmatrix.md),
  matching the weight that mean imputation already gives it in
  [`Gmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Gmatrix.md).
- [`alleleFrequency()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/alleleFrequency.md)
  no longer understates observed heterozygosity after imputation. Only
  genotype calls count towards it: the dosage must be a whole number and
  the cell must not be marked as imputed.
- A marker with no genotype calls at all is no longer filled in with a
  fabricated dosage of zero for every individual.
  [`imputeMarkers()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/imputeMarkers.md)
  drops such markers with a warning, controlled by the new
  `drop_uncalled` argument, and the relationship matrices drop them too
  instead of silently returning a matrix of `NaN`.
- [`BreedingExperiment()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment.md)
  matches `phenotypes` to individuals by row name when the row names are
  identifiers, instead of overwriting them by position. Row names that
  are neither absent nor matching are now an error. Duplicated genotype
  column names are rejected.
- A pedigree that gives the same individual different parents in
  different rows is now an error rather than being resolved by keeping
  whichever row came first; exact duplicate rows are dropped with a
  warning. Pedigree rows with a missing id are rejected.
- [`Hmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Hmatrix.md)
  warns when a genotyped individual does not appear in the pedigree,
  rather than dropping it silently.
- [`simulateBreeding()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/simulateBreeding.md)
  validates its arguments, and with `n_founder = 2` can no longer give
  an individual itself as both parents, which was a consequence of
  [`sample()`](https://rdrr.io/r/base/sample.html) treating a length-one
  vector as `1:n`. A population of a single individual is now allowed.

### New features

- [`imputeMarkers()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/imputeMarkers.md)
  records which cells it filled in a logical `imputed` assay, which
  travels with the object through subsetting and is what lets
  [`alleleFrequency()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/alleleFrequency.md)
  tell an imputed dosage from a call.

### Other changes

- [`simulateBreeding()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/simulateBreeding.md)
  now scales the simulated trait so that the genetic and residual
  variances are exactly `h2` and `1 - h2`, and `trueBV` is on that
  scale. `h2 = 0` and `h2 = 1` behave as they should. Simulated breeding
  values are therefore on a different scale than in earlier versions.

- The vignette opens with the motivation for the package and a
  comparison with `SNPRelate`, `GENESIS`, `GeneticsPed`, `Pedixplorer`,
  `FamAgg` and CRAN’s `AGHmatrix`, followed by an installation section.
  It now gives guidance on choosing `min_maf`, `min_call_rate`, `blend`
  and `tune`, and ends with a worked single-step BLUP: the mixed model
  equations solved with `H`, validated by predicting a withheld
  generation.

- Fixed a stray closing parenthesis in a vignette code block.

- Documented how phenotypes are matched to individuals, and what
  imputation costs in the marker summaries. \# BreedingExperiment 0.99.4

- Documented the return value on the class page, clearing the BiocCheck
  warning about a missing `alue` section.

- Wrapped the last over-long source line.

## BreedingExperiment 0.99.3

- Version bump only, to re-run the Bioconductor checks after a transient
  bioconductor.org outage (HTTP 504) aborted BiocCheck. No code changes.

## BreedingExperiment 0.99.2

- Split
  [`simulateBreeding()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/simulateBreeding.md)
  into helpers so no function exceeds the recommended length.
- Wrapped over-long source lines and added the `QualityControl`
  biocView.

## BreedingExperiment 0.99.1

- Removed a non-standard field from `DESCRIPTION`.
- Added citation metadata (`inst/CITATION` and `CITATION.cff`).

## BreedingExperiment 0.99.0

First submission to Bioconductor.

- The `BreedingExperiment` class, extending
  `RangedSummarizedExperiment`, holds marker genotypes, phenotypes, a
  pedigree and marker coordinates in one object.
- Relationship matrices:
  [`Amatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Amatrix.md)
  (pedigree),
  [`Gmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Gmatrix.md)
  (genomic),
  [`Dmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Dmatrix.md)
  (dominance) and
  [`Hmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Hmatrix.md)
  (single step, combining pedigree and genomic information so genotyped
  and ungenotyped individuals can be analysed together).
- Quality control with
  [`filterMarkers()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/filterMarkers.md)
  and
  [`imputeMarkers()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/imputeMarkers.md),
  plus
  [`alleleFrequency()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/alleleFrequency.md),
  [`inbreeding()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/inbreeding.md)
  and
  [`sortPedigree()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/sortPedigree.md).
- `demoBreeding`, a simulated five-generation population with
  ungenotyped ancestors, and
  [`simulateBreeding()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/simulateBreeding.md)
  for generating your own.
