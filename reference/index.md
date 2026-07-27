# Package index

## The container

A class that keeps genotypes, phenotypes, pedigree and marker
coordinates together, extending RangedSummarizedExperiment.

- [`show(`*`<BreedingExperiment>`*`)`](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment-class.md)
  : Breeding experiment container
- [`BreedingExperiment()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment.md)
  : Create a breeding experiment

## Relationship matrices

The matrices quantitative genetics runs on, from the pedigree, from the
markers, or from both at once.

- [`Amatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Amatrix.md)
  : Pedigree numerator relationship matrix
- [`Gmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Gmatrix.md)
  : Genomic relationship matrix
- [`Dmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Dmatrix.md)
  : Dominance genomic relationship matrix
- [`Hmatrix()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/Hmatrix.md)
  : Single-step relationship matrix
- [`inbreeding()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/inbreeding.md)
  : Inbreeding coefficients

## Accessors

- [`pedigree()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/pedigree.md)
  [`` `pedigree<-`() ``](https://mqfarooqi1.github.io/BreedingExperiment/reference/pedigree.md)
  : Pedigree of a breeding experiment
- [`genotypes()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/genotypes.md)
  [`` `genotypes<-`() ``](https://mqfarooqi1.github.io/BreedingExperiment/reference/genotypes.md)
  : Genotype dosages
- [`genotyped()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/genotyped.md)
  : Which individuals are genotyped
- [`alleleFrequency()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/alleleFrequency.md)
  : Allele frequencies and marker summaries

## Quality control

- [`filterMarkers()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/filterMarkers.md)
  : Filter markers on quality
- [`imputeMarkers()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/imputeMarkers.md)
  : Impute missing genotypes
- [`sortPedigree()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/sortPedigree.md)
  : Order a pedigree so that parents precede their offspring

## Data and simulation

- [`demoBreeding`](https://mqfarooqi1.github.io/BreedingExperiment/reference/demoBreeding.md)
  : A demonstration breeding population
- [`simulateBreeding()`](https://mqfarooqi1.github.io/BreedingExperiment/reference/simulateBreeding.md)
  : Simulate a small breeding population

## Package

- [`BreedingExperiment-package`](https://mqfarooqi1.github.io/BreedingExperiment/reference/BreedingExperiment-package.md)
  : BreedingExperiment: breeding data and relationship matrices
