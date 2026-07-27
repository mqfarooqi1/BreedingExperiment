# Allele frequencies and marker summaries

Computes, for every marker, the frequency of the counted allele, the
minor allele frequency, the call rate and the observed heterozygosity.

## Usage

``` r
alleleFrequency(x, ...)

# S4 method for class 'BreedingExperiment'
alleleFrequency(x, ...)
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
[S4Vectors::DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
with one row per marker and columns `freq`, `maf`, `callRate` and
`heterozygosity`.

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
```
