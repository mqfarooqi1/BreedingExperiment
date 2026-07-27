# BreedingExperiment: breeding data and relationship matrices

Keeps marker genotypes, phenotypes, pedigree and marker coordinates in
one object, and computes the relationship matrices that quantitative
genetics runs on.

## Details

The container,
[BreedingExperiment](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment-class.md),
extends
[SummarizedExperiment::RangedSummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html),
so markers are genomic ranges and the usual Bioconductor accessors
apply. The methods that matter are:

- [`Amatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Amatrix.md):

  pedigree numerator relationship matrix

- [`Gmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Gmatrix.md):

  genomic relationship matrix from markers

- [`Dmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Dmatrix.md):

  dominance relationship matrix

- [`Hmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Hmatrix.md):

  single-step matrix combining pedigree and genomic information, so
  genotyped and ungenotyped individuals can be analysed together

[`simulateBreeding()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/simulateBreeding.md)
produces a small pedigreed population for trying these out, and
[`filterMarkers()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/filterMarkers.md)
and
[`imputeMarkers()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/imputeMarkers.md)
cover routine quality control.

## Provenance

A substantial part of this package was written with the assistance of an
AI coding assistant (`Assisted-by: Claude, Anthropic`). The author
defined the scope, chose the statistical methods, and reviewed, tested
and validated the result, and is responsible for the correctness and
maintenance of the code.

Every method implemented here is an established, published technique in
quantitative genetics, cited in the documentation of the function that
implements it. The relationship matrices are checked in the test suite
against values that follow from theory rather than from a previous run:
a parent-offspring relationship of one half, a grandparent relationship
of one quarter, full sibs related by one half, the offspring of a
full-sib mating inbred at \\F = 0.25\\, a positive definite \\H\\, and
`Hmatrix(blend = 1)` reproducing \\A\_{22}\\ exactly.

## See also

Useful links:

- <https://github.com/mqfarooqi1/BreedingExperiment>

- <https://mqfarooqi1.github.io/BreedingExperiment/>

- Report bugs at
  <https://github.com/mqfarooqi1/BreedingExperiment/issues>

## Author

**Maintainer**: Muhammad Farooqi <mqfarooqi@gmail.com>
([ORCID](https://orcid.org/0000-0003-4918-9791))

Authors:

- Muhammad Farooqi <mqfarooqi@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-4918-9791))
