# Create a breeding experiment

Builds a
[BreedingExperiment](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment-class.md)
from a marker genotype matrix, with optional phenotypes, pedigree and
marker coordinates.

## Usage

``` r
BreedingExperiment(
  genotypes,
  phenotypes = NULL,
  pedigree = NULL,
  rowRanges = NULL,
  ...
)
```

## Arguments

- genotypes:

  A numeric matrix of allele dosages, markers in rows and individuals in
  columns, with column names identifying the individuals.

- phenotypes:

  Optional data frame or
  [S4Vectors::DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  of individual-level variables, with one row per column of `genotypes`.

- pedigree:

  Optional data frame with columns `id`, `sire` and `dam`. Unknown
  parents may be `NA` or `0`.

- rowRanges:

  Optional
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  giving the position of each marker, of the same length as
  `nrow(genotypes)`.

- ...:

  Further assays, passed to
  [`SummarizedExperiment::SummarizedExperiment()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html).

## Value

A
[BreedingExperiment](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment-class.md)
object.

## Details

Genotypes are allele dosages: usually the integers 0, 1 and 2, with `NA`
for missing calls. Fractional values in between are accepted, so
expected dosages from imputation can be stored directly. Markers are
rows and individuals are columns, matching the Bioconductor convention
that features are rows.

If `rowRanges` is not supplied, markers are given placeholder
coordinates on a single sequence so that the object is still a valid
ranged experiment; pass real coordinates whenever you have them, because
that is what allows a data set to be subset by genomic region.

The pedigree may name ancestors that were never genotyped. Those extra
individuals are exactly what
[`Hmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Hmatrix.md)
needs in order to combine pedigree and genomic information.

## See also

[`Gmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Gmatrix.md),
[`Amatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Amatrix.md),
[`Hmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Hmatrix.md),
[`simulateBreeding()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/simulateBreeding.md)

## Examples

``` r
set.seed(1)
geno <- matrix(rbinom(40 * 12, 2, 0.3), nrow = 40, ncol = 12,
               dimnames = list(paste0("snp", 1:40), paste0("ind", 1:12)))
pheno <- data.frame(yield = rnorm(12))
be <- BreedingExperiment(geno, phenotypes = pheno)
be
#> class: BreedingExperiment
#> markers: 40  individuals genotyped: 12 
#> assays(1): genotype
#> phenotypes(1): yield
#> sequences(1): unknown
#> pedigree: none
dim(be)
#> [1] 40 12
```
