# Genotype dosages

Gets or sets the matrix of allele dosages, a shorthand for the
`genotype` assay. Markers are rows and individuals are columns.

## Usage

``` r
genotypes(x, ...)

genotypes(x) <- value

# S4 method for class 'BreedingExperiment'
genotypes(x, ...)

# S4 method for class 'BreedingExperiment'
genotypes(x) <- value
```

## Arguments

- x:

  A
  [BreedingExperiment](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment-class.md)
  object.

- ...:

  Unused.

- value:

  A numeric matrix of dosages coded 0, 1 or 2.

## Value

`genotypes()` returns a numeric matrix; the replacement form returns the
updated object.

## Examples

``` r
set.seed(1)
be <- simulateBreeding(n_ind = 6, n_marker = 10)
genotypes(be)[1:5, 1:4]
#>         ind001 ind002 ind003 ind004
#> snp0001      0      2      1      1
#> snp0002      0      0      1      1
#> snp0003      1      1      1      1
#> snp0004      2      2      2      2
#> snp0005      1      0      1      0
```
