# Impute missing genotypes

Fills missing dosages with twice the allele frequency of the marker, the
simple mean imputation that relationship matrices assume. It is adequate
for the low rates of missingness left after quality control; for
anything more, impute with dedicated software before building the
object.

## Usage

``` r
imputeMarkers(x, ...)

# S4 method for class 'BreedingExperiment'
imputeMarkers(x, ...)
```

## Arguments

- x:

  A
  [BreedingExperiment](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment-class.md)
  object.

- ...:

  Unused.

## Value

A
[BreedingExperiment](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment-class.md)
object with no missing genotypes.

## See also

[`filterMarkers()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/filterMarkers.md)

## Examples

``` r
set.seed(1)
be <- simulateBreeding(n_ind = 20, n_marker = 50, missing = 0.1)
anyNA(genotypes(be))
#> [1] TRUE
anyNA(genotypes(imputeMarkers(be)))
#> [1] FALSE
```
