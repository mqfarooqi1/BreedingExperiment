# Allele frequencies and marker summaries

Computes, for every marker, the frequency of the counted allele, the
minor allele frequency, the call rate and the observed heterozygosity.

## Usage

``` r
alleleFrequency(x, ...)

# S4 method for class 'BreedingExperiment'
alleleFrequency(x, tol = 1e-08, ...)
```

## Arguments

- x:

  A
  [BreedingExperiment](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment-class.md)
  object.

- ...:

  Unused.

- tol:

  How far a dosage may sit from a whole number and still count as a
  genotype call when `heterozygosity` is computed.

## Value

A
[S4Vectors::DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
with one row per marker and columns `freq`, `maf`, `callRate` and
`heterozygosity`.

## Details

`freq` and `maf` use every dosage that is present, so an expected dosage
from imputation contributes to them on the same footing as a call.
`callRate` is the proportion of dosages that are not `NA`; after
[`imputeMarkers()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/imputeMarkers.md)
that is 1 for every marker, so take call rates before imputing if you
need them.

`heterozygosity` is different, because a heterozygote is a genotype
rather than a quantity. A cell counts towards it only if it is a
genotype call: its dosage must be a whole number, and it must not be
marked in the `imputed` assay that
[`imputeMarkers()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/imputeMarkers.md)
writes. Everything else is left out of both the numerator and the
denominator. Without those two checks an imputed cell would be counted
as a non-heterozygote and heterozygosity would appear to fall every time
a data set was imputed. A marker with no genotype call at all gets `NA`.

Dosages imputed by other software arrive without the `imputed` assay, so
only the whole-number check applies to them; pass such dosages in as
fractional values, or carry your own mask as an extra assay, if observed
heterozygosity matters to you.

## See also

[`filterMarkers()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/filterMarkers.md),
[`imputeMarkers()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/imputeMarkers.md)

## Examples

``` r
set.seed(1)
be <- simulateBreeding(n_ind = 30, n_marker = 50)
head(alleleFrequency(be))
#> DataFrame with 6 rows and 4 columns
#>              freq       maf  callRate heterozygosity
#>         <numeric> <numeric> <numeric>      <numeric>
#> snp0001  0.450000  0.450000         1       0.633333
#> snp0002  0.333333  0.333333         1       0.466667
#> snp0003  0.500000  0.500000         1       0.600000
#> snp0004  0.816667  0.183333         1       0.233333
#> snp0005  0.333333  0.333333         1       0.333333
#> snp0006  0.800000  0.200000         1       0.266667

# imputation does not change the observed heterozygosity of a marker,
# because imputed cells are not counted as genotype calls
bem <- simulateBreeding(n_ind = 30, n_marker = 50, missing = 0.2)
all.equal(alleleFrequency(bem)$heterozygosity,
          alleleFrequency(imputeMarkers(bem))$heterozygosity)
#> [1] TRUE
```
