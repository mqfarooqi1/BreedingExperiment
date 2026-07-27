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

  Number of unrelated founders.

- genotyped:

  Number of individuals with genotypes, taken from the most recent ones.
  Defaults to all of them.

- n_qtl:

  Number of markers with an effect on the trait.

- h2:

  Narrow-sense heritability of the simulated trait.

- missing:

  Proportion of genotype calls set to `NA`.

## Value

A
[BreedingExperiment](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment-class.md)
object with a `phenotype` column and the true breeding value in
`colData`.

## Details

Set `genotyped` below `n_ind` to leave some individuals ungenotyped,
which is the situation
[`Hmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Hmatrix.md)
is designed for.

The function draws random numbers but does not set the seed itself. Call
[`set.seed()`](https://rdrr.io/r/base/Random.html) beforehand if you
need the same population twice.

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
#> ind001  0.04962359  0.8492997
#> ind002  2.09499329 -0.4066208
#> ind003 -0.69175387 -0.4045620
#> ind004 -1.27278074 -0.1810128
#> ind005 -1.83665981 -0.1741050
#> ind006 -1.12827745 -0.4129922
```
