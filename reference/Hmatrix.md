# Single-step relationship matrix

Builds the matrix \\H\\ that combines pedigree and genomic information,
so that genotyped and ungenotyped individuals can be analysed in a
single model. This is the relationship structure behind single-step
genomic prediction (ssGBLUP).

## Usage

``` r
Hmatrix(x, ...)

# S4 method for class 'BreedingExperiment'
Hmatrix(x, blend = 0.05, tune = TRUE, min_maf = 0, impute = TRUE, ...)
```

## Arguments

- x:

  A
  [BreedingExperiment](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment-class.md)
  object with both a pedigree and genotypes.

- ...:

  Unused.

- blend:

  Weight \\w\\ given to \\A\_{22}\\ when blending, between 0 and 1. The
  default, 0.05, is a common choice.

- tune:

  Rescale \\G\\ to the scale of \\A\_{22}\\ before blending.

- min_maf, impute:

  Passed to
  [`Gmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Gmatrix.md).

## Value

A symmetric numeric matrix over every individual in the pedigree,
ungenotyped individuals first, with an attribute `genotyped` giving the
identifiers of the genotyped group.

## Details

Real breeding programmes genotype only part of the population. Using
only \\G\\ throws away every ungenotyped relative; using only \\A\\
throws away the markers. \\H\\ keeps both: relationships among genotyped
individuals come from the markers, relationships among ungenotyped
individuals are those from the pedigree corrected by what the markers
reveal about their genotyped relatives, and the two groups are connected
through the pedigree.

Writing 1 for ungenotyped and 2 for genotyped individuals, \$\$H\_{22} =
G^{\*}\$\$ \$\$H\_{12} = A\_{12} A\_{22}^{-1} G^{\*}\$\$ \$\$H\_{11} =
A\_{11} + A\_{12} A\_{22}^{-1} (G^{\*} - A\_{22}) A\_{22}^{-1}
A\_{21}\$\$

Because \\G\\ and \\A\_{22}\\ are estimated on different scales, \\G\\
is first made compatible with \\A\_{22}\\ in two standard steps.
*Tuning* rescales \\G\\ so that its mean diagonal and mean off-diagonal
match those of \\A\_{22}\\. *Blending* then mixes in a small share of
\\A\_{22}\\, \\G^{\*} = (1-w)G + wA\_{22}\\, which also makes the result
invertible when markers are fewer than individuals.

## References

Legarra, A., Aguilar, I. & Misztal, I. (2009) "A relationship matrix
including full pedigree and genomic information." Journal of Dairy
Science 92, 4656-4663.
[doi:10.3168/jds.2009-2061](https://doi.org/10.3168/jds.2009-2061)

Christensen, O. F. & Lund, M. S. (2010) "Genomic prediction when some
animals are not genotyped." Genetics Selection Evolution 42, 2.
[doi:10.1186/1297-9686-42-2](https://doi.org/10.1186/1297-9686-42-2)

Aguilar, I., Misztal, I., Johnson, D. L., Legarra, A., Tsuruta, S. &
Lawlor, T. J. (2010) "Hot topic: A unified approach to utilize
phenotypic, full pedigree, and genomic information for genetic
evaluation of Holstein final score." Journal of Dairy Science 93,
743-752.
[doi:10.3168/jds.2009-2730](https://doi.org/10.3168/jds.2009-2730)

## See also

[`Amatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Amatrix.md),
[`Gmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Gmatrix.md)

## Examples

``` r
# only some of the individuals are genotyped
set.seed(1)
be <- simulateBreeding(n_ind = 30, n_marker = 300, genotyped = 20)
H <- Hmatrix(be)
dim(H)
#> [1] 30 30
attr(H, "genotyped")
#>  [1] "ind011" "ind013" "ind019" "ind024" "ind027" "ind030" "ind012" "ind014"
#>  [9] "ind020" "ind022" "ind026" "ind015" "ind016" "ind017" "ind018" "ind021"
#> [17] "ind023" "ind025" "ind028" "ind029"
```
