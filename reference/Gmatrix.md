# Genomic relationship matrix

Builds the additive genomic relationship matrix \\G\\ from marker
dosages using VanRaden's first method: with \\W\\ the marker matrix
centred by twice the allele frequency, \\G = WW' / 2\sum p(1-p)\\.

## Usage

``` r
Gmatrix(x, ...)

# S4 method for class 'BreedingExperiment'
Gmatrix(x, min_maf = 0, impute = TRUE, ...)
```

## Arguments

- x:

  A
  [BreedingExperiment](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment-class.md)
  object.

- ...:

  Unused.

- min_maf:

  Markers with a minor allele frequency below this are dropped.

- impute:

  Replace missing dosages by twice the allele frequency before
  computing. If `FALSE`, missing values raise an error.

## Value

A symmetric numeric matrix over the genotyped individuals.

## Details

\\G\\ plays the same role as the pedigree matrix \\A\\ but is estimated
from markers, so it captures Mendelian sampling that a pedigree cannot:
full sibs have an expected pedigree relationship of exactly 0.5, whereas
their genomic relationship varies around it.

## References

VanRaden, P. M. (2008) "Efficient methods to compute genomic
predictions." Journal of Dairy Science 91, 4414-4423.
[doi:10.3168/jds.2007-0980](https://doi.org/10.3168/jds.2007-0980)

## See also

[`Amatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Amatrix.md),
[`Dmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Dmatrix.md),
[`Hmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Hmatrix.md)

## Examples

``` r
set.seed(1)
be <- simulateBreeding(n_ind = 20, n_marker = 200)
G <- Gmatrix(be)
round(G[1:4, 1:4], 3)
#>        ind001 ind002 ind003 ind004
#> ind001  0.950 -0.108 -0.107 -0.061
#> ind002 -0.108  1.049 -0.025  0.008
#> ind003 -0.107 -0.025  0.930 -0.125
#> ind004 -0.061  0.008 -0.125  1.198
```
