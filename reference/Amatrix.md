# Pedigree numerator relationship matrix

Builds the additive relationship matrix \\A\\ from the pedigree using
the tabular method. Element \\A\_{ij}\\ is twice the coefficient of
kinship between individuals \\i\\ and \\j\\, so a diagonal entry is
\\1 + F_i\\ where \\F_i\\ is the inbreeding coefficient.

## Usage

``` r
Amatrix(x, ...)

# S4 method for class 'BreedingExperiment'
Amatrix(x, ids = NULL, ...)
```

## Arguments

- x:

  A
  [BreedingExperiment](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment-class.md)
  object.

- ...:

  Unused.

- ids:

  Optional character vector of identifiers to return, in that order.
  Defaults to every individual in the pedigree.

## Value

A symmetric numeric matrix with dimnames taken from the pedigree.

## Details

The matrix is built recursively over a pedigree sorted so that parents
come first: an individual's relationship to an earlier individual is the
average of its parents' relationships to that individual, and its
diagonal is \\1 + \tfrac{1}{2}A\_{sire,dam}\\.

## References

Henderson, C. R. (1976) "A simple method for computing the inverse of a
numerator relationship matrix used in prediction of breeding values."
Biometrics 32, 69-83.
[doi:10.2307/2529339](https://doi.org/10.2307/2529339)

## See also

[`Gmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Gmatrix.md),
[`Hmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Hmatrix.md),
[`inbreeding()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/inbreeding.md)

## Examples

``` r
set.seed(1)
be <- simulateBreeding(n_ind = 12, n_marker = 30)
A <- Amatrix(be)
round(A[1:5, 1:5], 3)
#>        ind001 ind002 ind003 ind004 ind005
#> ind001      1      0      0      0      0
#> ind002      0      1      0      0      0
#> ind003      0      0      1      0      0
#> ind004      0      0      0      1      0
#> ind005      0      0      0      0      1
```
