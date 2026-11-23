# Simulate a small breeding population

Generates a pedigreed population with markers and a heritable trait, for
examples, tests and teaching. Founders are unrelated; later generations
are produced by mating sampled parents, and each offspring inherits one
allele from each parent at every marker, so the simulated genotypes and
the pedigree agree with each other.

## Usage

``` r
simulateBreeding(
  n_ind = 50,
  n_marker = 200,
  n_founder = 10,
  genotyped = NULL,
  n_qtl = 20,
  h2 = 0.4,
  missing = 0
)
```

## Arguments

- n_ind:

  Number of individuals in total.

- n_marker:

  Number of markers.

- n_founder:

  Number of unrelated founders. Capped at `n_ind`; at least two are
  needed if there are to be any offspring.

- genotyped:

  Number of individuals with genotypes, taken from the most recent ones.
  Defaults to all of them.

- n_qtl:

  Number of markers with an effect on the trait.

- h2:

  Narrow-sense heritability of the simulated trait, between 0 and 1.

- missing:

  Proportion of genotype calls set to `NA`, between 0 and 1.

## Value

A
[BreedingExperiment](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment-class.md)
object with a `phenotype` column and the true breeding value in
`colData`.

## Details

Set `genotyped` below `n_ind` to leave some individuals ungenotyped,
which is the situation
[`Hmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Hmatrix.md)
is designed for. The genotyped individuals are the most recent ones, so
the ungenotyped are the earlier generations, as in a programme that
began genotyping part way through.

Note that `n_founder` defaults to 10 and is capped at `n_ind`, so a call
such as `simulateBreeding(n_ind = 8)` returns eight unrelated founders
and no pedigree structure at all. Ask for more individuals than founders
to get relatives.

The trait is scaled so that the genetic variance is `h2` and the
residual variance is `1 - h2`, and `trueBV` in the column data is the
genetic component on that scale.

The function draws random numbers but does not set the seed itself. Call
[`set.seed()`](https://rdrr.io/r/base/Random.html) beforehand if you
need the same population twice.

## See also

[`BreedingExperiment()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment.md),
[`Hmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Hmatrix.md)

## Examples

``` r
set.seed(1)
be <- simulateBreeding(n_ind = 25, n_marker = 100)
be
#> class: BreedingExperiment
#> markers: 100  individuals genotyped: 25 
#> assays(1): genotype
#> phenotypes(2): phenotype trueBV
#> sequences(5): chr1 chr2 chr3 chr4 chr5
#> pedigree: 25 individuals (25 genotyped, 0 not; 10 founders)
head(as.data.frame(colData(be)))
#>          phenotype     trueBV
#> ind001  0.03138471  0.5371443
#> ind002  1.32499010 -0.2571695
#> ind003 -0.43750356 -0.2558675
#> ind004 -0.80497722 -0.1144826
#> ind005 -1.16160565 -0.1101137
#> ind006 -0.71358531 -0.2611992

# a population where only the recent individuals were genotyped
part <- simulateBreeding(n_ind = 30, n_marker = 100, n_founder = 6,
                         genotyped = 20)
c(pedigree = nrow(pedigree(part)), genotyped = ncol(part))
#>  pedigree genotyped 
#>        30        20 
```
