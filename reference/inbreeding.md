# Inbreeding coefficients

Returns the inbreeding coefficient of each individual, obtained from the
diagonal of the pedigree relationship matrix as \\F_i = A\_{ii} - 1\\.

## Usage

``` r
inbreeding(x, ...)

# S4 method for class 'BreedingExperiment'
inbreeding(x, ...)
```

## Arguments

- x:

  A
  [BreedingExperiment](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment-class.md)
  object.

- ...:

  Unused.

## Value

A named numeric vector of inbreeding coefficients.

## See also

[`Amatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Amatrix.md)

## Examples

``` r
set.seed(1)
be <- simulateBreeding(n_ind = 12, n_marker = 30)
round(inbreeding(be), 4)
#> ind001 ind002 ind003 ind004 ind005 ind006 ind007 ind008 ind009 ind010 ind011 
#>      0      0      0      0      0      0      0      0      0      0      0 
#> ind012 
#>      0 
```
