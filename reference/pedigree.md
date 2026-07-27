# Pedigree of a breeding experiment

Gets or sets the pedigree, a table with columns `id`, `sire` and `dam`.
The pedigree may include ancestors that were never genotyped; those
individuals are what allow
[`Hmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Hmatrix.md)
to combine pedigree and genomic information.

## Usage

``` r
pedigree(x, ...)

pedigree(x) <- value

# S4 method for class 'BreedingExperiment'
pedigree(x, ...)

# S4 method for class 'BreedingExperiment'
pedigree(x) <- value
```

## Arguments

- x:

  A
  [BreedingExperiment](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment-class.md)
  object.

- ...:

  Unused.

- value:

  A data frame with columns `id`, `sire` and `dam`.

## Value

`pedigree()` returns a
[S4Vectors::DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html);
the replacement form returns the updated object.

## Examples

``` r
set.seed(1)
be <- simulateBreeding(n_ind = 8, n_marker = 20)
head(pedigree(be))
#> DataFrame with 6 rows and 3 columns
#>            id        sire         dam
#>   <character> <character> <character>
#> 1      ind001          NA          NA
#> 2      ind002          NA          NA
#> 3      ind003          NA          NA
#> 4      ind004          NA          NA
#> 5      ind005          NA          NA
#> 6      ind006          NA          NA
```
