# Dominance genomic relationship matrix

Builds the dominance relationship matrix from marker dosages using the
genotypic parameterisation of Vitezica, Varona and Legarra (2013).
Heterozygotes and the two homozygotes are coded so that the resulting
matrix captures dominance deviations, orthogonal to the additive
relationships returned by
[`Gmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Gmatrix.md).

## Usage

``` r
Dmatrix(x, ...)

# S4 method for class 'BreedingExperiment'
Dmatrix(x, min_maf = 0, impute = TRUE, ...)
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

Fitting additive and dominance effects together lets a model separate
the part of the genetic variance that is inherited predictably from the
part that is not, which matters for crossbreeding and for predicting the
performance of specific crosses.

## References

Vitezica, Z. G., Varona, L. & Legarra, A. (2013) "On the additive and
dominant variance and covariance of individuals within the genomic
selection scope." Genetics 195, 1223-1230.
[doi:10.1534/genetics.113.155176](https://doi.org/10.1534/genetics.113.155176)

## See also

[`Gmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Gmatrix.md)

## Examples

``` r
set.seed(1)
be <- simulateBreeding(n_ind = 20, n_marker = 200)
D <- Dmatrix(be)
round(D[1:4, 1:4], 3)
#>        ind001 ind002 ind003 ind004
#> ind001  0.865 -0.010  0.133 -0.033
#> ind002 -0.010  1.085  0.019  0.055
#> ind003  0.133  0.019  1.021 -0.100
#> ind004 -0.033  0.055 -0.100  1.136
```
