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

  Allow missing dosages, giving them a dominance contribution of zero.
  If `FALSE`, missing values raise an error.

## Value

A symmetric numeric matrix over the genotyped individuals.

## Details

Fitting additive and dominance effects together lets a model separate
the part of the genetic variance that is inherited predictably from the
part that is not, which matters for crossbreeding and for predicting the
performance of specific crosses.

## Dosages that are not whole numbers

Vitezica and others code the genotypes AA, Aa and aa as \\-2q^2\\,
\\2pq\\ and \\-2p^2\\. Those three values lie on a single quadratic in
the dosage \\g\\, \$\$-g^2 + (2p + 1)g - 2p^2,\$\$ which is what this
function evaluates. At dosages 2, 1 and 0 it returns the three published
values exactly, and an expected dosage from imputation falls between
them rather than being forced into one genotype class.

A missing call is given a dominance contribution of zero, which is the
weight that mean imputation already gives it in
[`Gmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Gmatrix.md).
A cell that
[`imputeMarkers()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/imputeMarkers.md)
has filled is treated the same way, because an expected dosage is not a
genotype: at a dosage of \\2p\\ the quadratic returns exactly \\2pq\\,
so counting it would record every imputed cell as a heterozygote. The
result is therefore the same whether or not
[`imputeMarkers()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/imputeMarkers.md)
has been run, provided the object carries the `imputed` assay that
records which cells were filled.

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
