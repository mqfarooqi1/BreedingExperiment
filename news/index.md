# Changelog

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
