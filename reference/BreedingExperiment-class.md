# Breeding experiment container

A `BreedingExperiment` keeps everything a quantitative-genetic analysis
needs in one object: marker genotypes, the phenotypes and design
variables of the individuals, the pedigree, and the position of every
marker on the genome.

## Usage

``` r
# S4 method for class 'BreedingExperiment'
show(object)
```

## Arguments

- object:

  A BreedingExperiment object.

## Value

A class definition, so nothing is returned by the class itself. The
constructor
[`BreedingExperiment()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment.md)
returns a `BreedingExperiment` object, and the `show` method is called
for the summary it prints, returning `object` invisibly.

## Details

The class extends
[SummarizedExperiment::RangedSummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html),
so markers are rows and individuals are columns, and the familiar
accessors work as usual:
[`assay()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
returns the genotype matrix, `rowRanges()` the marker coordinates as a
[GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html),
and
[`colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
the phenotypes. Because marker positions are genomic ranges, a breeding
data set can be subset by region, overlapped with annotation, and
otherwise handled with the standard Bioconductor vocabulary.

One slot is added to the parent class:

- `pedigree`:

  A
  [S4Vectors::DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  with columns `id`, `sire` and `dam`, giving the parents of each
  individual. Unknown parents are `NA`. It may contain ancestors that
  were never genotyped, which is what makes single-step analysis
  possible.

Genotypes are stored as allele dosages: 0, 1 or 2 copies of the counted
allele, with `NA` for missing calls.

## Construction

Use
[`BreedingExperiment()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment.md);
see its help page for examples.

## Display

The `show` method prints the size of the object, the assays it holds,
the phenotype columns, the sequences the markers lie on, and a summary
of the pedigree including how many individuals are genotyped.

## See also

[`BreedingExperiment()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment.md),
[`Amatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Amatrix.md),
[`Gmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Gmatrix.md),
[`Hmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Hmatrix.md)

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
```
