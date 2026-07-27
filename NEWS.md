# BreedingExperiment 0.99.3

* Version bump only, to re-run the Bioconductor checks after a transient
  bioconductor.org outage (HTTP 504) aborted BiocCheck. No code changes.

# BreedingExperiment 0.99.2

* Split `simulateBreeding()` into helpers so no function exceeds the recommended length.
* Wrapped over-long source lines and added the `QualityControl` biocView.

# BreedingExperiment 0.99.1

* Removed a non-standard field from `DESCRIPTION`.
* Added citation metadata (`inst/CITATION` and `CITATION.cff`).

# BreedingExperiment 0.99.0

First submission to Bioconductor.

* The `BreedingExperiment` class, extending `RangedSummarizedExperiment`, holds
  marker genotypes, phenotypes, a pedigree and marker coordinates in one object.
* Relationship matrices: `Amatrix()` (pedigree), `Gmatrix()` (genomic),
  `Dmatrix()` (dominance) and `Hmatrix()` (single step, combining pedigree and
  genomic information so genotyped and ungenotyped individuals can be analysed
  together).
* Quality control with `filterMarkers()` and `imputeMarkers()`, plus
  `alleleFrequency()`, `inbreeding()` and `sortPedigree()`.
* `demoBreeding`, a simulated five-generation population with ungenotyped
  ancestors, and `simulateBreeding()` for generating your own.
