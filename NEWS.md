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
