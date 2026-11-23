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
imputeMarkers(x, drop_uncalled = TRUE, ...)
```

## Arguments

- x:

  A
  [BreedingExperiment](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment-class.md)
  object.

- ...:

  Unused.

- drop_uncalled:

  Drop markers with no genotype calls at all, with a warning. If
  `FALSE`, their presence is an error.

## Value

A
[BreedingExperiment](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment-class.md)
object with no missing genotypes and an `imputed` assay marking the
cells that were filled in.

## Details

Mean imputation leaves an expected dosage rather than a whole number of
copies, which the container accepts. Note what that costs: an imputed
cell is no longer a genotype call, and its value alone no longer says
so, since a marker whose mean dosage is exactly 1 imputes to exactly 1.
The cells that were filled in are therefore recorded in a logical assay
named `imputed`, reachable with `assay(x, "imputed")` and carried
through subsetting with everything else.
[`alleleFrequency()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/alleleFrequency.md)
consults it, so observed heterozygosity is unchanged by imputing.
`callRate`, on the other hand, counts dosages that are present and so
becomes 1 for every marker; record call rates before imputing if you
need them.

A marker with no calls at all has no frequency to impute from. Such
markers are dropped with a warning rather than being filled in, because
filling them would invent a genotype for every individual at that
marker.
[`filterMarkers()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/filterMarkers.md)
removes them too, and running it first is the better habit.

## See also

[`filterMarkers()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/filterMarkers.md),
[`alleleFrequency()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/alleleFrequency.md)

## Examples

``` r
set.seed(1)
be <- simulateBreeding(n_ind = 20, n_marker = 50, missing = 0.1)
anyNA(genotypes(be))
#> [1] TRUE
anyNA(genotypes(imputeMarkers(be)))
#> [1] FALSE
```
