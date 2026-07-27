# Filter markers on quality

Removes markers that are uninformative or poorly genotyped: those below
a minor allele frequency threshold, those called in too few individuals,
and optionally those that are monomorphic.

## Usage

``` r
filterMarkers(x, ...)

# S4 method for class 'BreedingExperiment'
filterMarkers(
  x,
  min_maf = 0.01,
  min_call_rate = 0.9,
  drop_monomorphic = TRUE,
  ...
)
```

## Arguments

- x:

  A
  [BreedingExperiment](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment-class.md)
  object.

- ...:

  Unused.

- min_maf:

  Minimum minor allele frequency to keep.

- min_call_rate:

  Minimum proportion of individuals with a call.

- drop_monomorphic:

  Drop markers with no variation.

## Value

A
[BreedingExperiment](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment-class.md)
object with fewer rows.

## See also

[`imputeMarkers()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/imputeMarkers.md),
[`alleleFrequency()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/alleleFrequency.md)

## Examples

``` r
set.seed(1)
be <- simulateBreeding(n_ind = 30, n_marker = 100, missing = 0.05)
nrow(be)
#> [1] 100
nrow(filterMarkers(be, min_maf = 0.05, min_call_rate = 0.9))
#> [1] 90
```
