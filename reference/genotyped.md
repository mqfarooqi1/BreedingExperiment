# Which individuals are genotyped

Returns the identifiers that appear in the genotype matrix. In a
single-step analysis the pedigree usually names more individuals than
this, and the difference is what
[`Hmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Hmatrix.md)
reconciles.

## Usage

``` r
genotyped(x, ...)

# S4 method for class 'BreedingExperiment'
genotyped(x, ...)
```

## Arguments

- x:

  A
  [BreedingExperiment](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment-class.md)
  object.

- ...:

  Unused.

## Value

A character vector of individual identifiers.

## Examples

``` r
set.seed(1)
be <- simulateBreeding(n_ind = 8, n_marker = 20)
genotyped(be)
#> [1] "ind001" "ind002" "ind003" "ind004" "ind005" "ind006" "ind007" "ind008"
```
