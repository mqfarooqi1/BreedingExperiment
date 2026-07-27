# A demonstration breeding population

A small, five-generation breeding population supplied so that every
function in the package can be tried on data that behaves like a real
programme.

## Usage

``` r
data(demoBreeding)
```

## Format

A
[BreedingExperiment](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment-class.md)
with 900 markers and 120 genotyped individuals, and a pedigree of 180.
`colData` holds:

- generation:

  Generation number, 3 to 5 for genotyped individuals.

- sex:

  `"F"` or `"M"`.

- yield:

  A moderately heritable trait (\\h^2 \approx 0.35\\), recorded on
  females only.

- stature:

  A more heritable trait (\\h^2 \approx 0.60\\), recorded on everyone.

- trueBV_yield, trueBV_stature:

  The simulated true breeding values.

## Source

Simulated by `data-raw/make_demo.R` in the package sources.

## Details

The data are **simulated**, not observed. They are built to have the
awkward features that real breeding data has and that tidy textbook
examples lack:

- **Only part of the population is genotyped.** The two earliest
  generations were born before genotyping began, so 60 of the 180
  individuals in the pedigree have no markers. This is the situation
  [`Hmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Hmatrix.md)
  exists for.

- **Unequal family sizes**, because a few sires are used heavily.

- **A trait recorded on one sex only.** `yield` is missing for males, as
  a milk trait would be.

- **Missing genotype calls**, at about one per cent, as any chip has.

- **Markers on ten chromosomes** with irregular spacing, so subsetting
  by genomic region is meaningful.

Because the data are simulated, the true breeding values are known and
are kept in `colData`, which makes the object useful for checking that a
prediction method recovers what it should.

## See also

[`simulateBreeding()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/simulateBreeding.md)
to generate your own,
[`Hmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Hmatrix.md)
for the analysis this data set is shaped for.

## Examples

``` r
data(demoBreeding)
demoBreeding
#> class: BreedingExperiment
#> markers: 900  individuals genotyped: 120 
#> assays(1): genotype
#> phenotypes(6): generation sex yield stature trueBV_yield trueBV_stature
#> sequences(10): chr1 chr2 chr3 chr4 chr5 chr6
#> pedigree: 180 individuals (120 genotyped, 60 not; 30 founders)

# 60 of the 180 individuals in the pedigree were never genotyped
nrow(pedigree(demoBreeding)) - ncol(demoBreeding)
#> [1] 60

# the single-step matrix covers all of them
H <- Hmatrix(demoBreeding)
dim(H)
#> [1] 180 180
```
