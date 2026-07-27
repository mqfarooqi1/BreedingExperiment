# Order a pedigree so that parents precede their offspring

The tabular method for building a relationship matrix requires that
every individual appears after both of its parents. This performs that
topological sort and reports a cycle if the pedigree contains one.

## Usage

``` r
sortPedigree(ped)
```

## Arguments

- ped:

  A data frame or
  [S4Vectors::DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  with columns `id`, `sire` and `dam`.

## Value

A data frame with the same columns, reordered.

## Examples

``` r
ped <- data.frame(id = c("c", "a", "b"),
                  sire = c("a", NA, NA),
                  dam = c("b", NA, NA))
sortPedigree(ped)
#>   id sire  dam
#> 2  a <NA> <NA>
#> 3  b <NA> <NA>
#> 1  c    a    b
```
