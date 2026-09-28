# Comparing SRNs

This repository contains the R implementation for the paper **"Stochastic ordering tools for continuous-time Markov chains and applications to reaction network models"** by Daniele Cappelletti, Giulio Cuniberti, and Paola Siri (https://doi.org/10.48550/arXiv.2604.00756). The code consists of two main scripts:
- `Examples.R`
- `Functions.R`

## Installation and usage

The code requires R and the packages `rcdd` and `gmp`, which provide exact rational and integer arithmetic. Install the packages once in R:

```r
install.packages(c("rcdd", "gmp"))
```

With the repository folder as your working directory, load the functions:

```r
source("Functions.R")
```

After sourcing `Functions.R`, you can use its two main functions:

#### `checkOrdering(source.complexes, product.complexes, preorder.matrix)`
checks whether Theorem 4.1 can be applied to the reaction network defined by the first two arguments, given the preorder specified by the third argument.

#### `findOrderings(source.complexes, product.complexes)`
searches the simple preorders considered in the paper, comparing each species count by equality, either inequality, or no constraint. It returns the associated rate inequalities for which Theorem 4.1 applies, identifying comparisons obtained from one another by reversing all inequalities.

The conditions are sufficient: an unsuccessful check or search does not prove that the processes cannot be ordered.

## Input format

- `source.complexes` and `product.complexes` are `n × d` matrices, where:
  - `n` = number of reactions
  - `d` = number of species
  - row `r` in `source.complexes` represents the **source complex** of reaction `r`
  - row `r` in `product.complexes` represents the **product complex** of the same reaction `r`
  - column `s` refers to the same species `s`, in both matrices

- `preorder.matrix` is the matrix *M* from Theorem 4.1

All matrix entries must be integers; source and product complexes must have nonnegative entries. Inputs must have absolute value less than `2^53`; computed integers outside the supported exact range are rejected. The kinetic and initial-state assumptions of Theorem 4.1 must also hold when applying its conclusion.

## Output description

**`checkOrdering(source.complexes, product.complexes, preorder.matrix)`** returns a list of 4 objects.
- `result`: 1 or 0, depending on whether the hypotheses of Theorem 4.1 are satisfied or not
- `rate.inequalities`: vector whose component `r` corresponds to reaction `r` and can take the values
  - `"="` if the rate constants are equal in *X* and *Y*
  - `"+"` if the rate constant of *Y* is greater than or equal to that of *X*
  - `"-"` if the rate constant of *Y* is less than or equal to that of *X*
  - `"?"` if the rate constants of *X* and *Y* are not compared
- `species.inequalities`: vector whose component `s` corresponds to species `s` and can take the values
  - `"="` if the molecular counts are equal in *X* and *Y*
  - `"+"` if the molecular count of *Y* is greater than or equal to that of *X*
  - `"-"` if the molecular count of *Y* is less than or equal to that of *X*
  - `"?"` if the molecular counts of *X* and *Y* are not compared
- `first.fail`: when `result = 0`, it specifies the first hypothesis of Theorem 4.1 that failed to be satisfied

The reported inequalities should be used only when `result = 1`.

**`findOrderings(source.complexes, product.complexes)`** returns a list of 2 objects.
- `preorders`: list of all pairs of `rate.inequalities` and `species.inequalities` (as described for the previous function) that define an admissible preordering structure
  and for which at least one component of `species.inequalities` is either `"+"` or `"-"`
- `equivalences`: list of all pairs of `rate.inequalities` and `species.inequalities` (as described for the previous function) that define the remaining admissible preordering (equivalence) structures

The two trivial comparisons (identical processes, or no species counts compared) are omitted from these lists.

## Examples

All examples presented in the paper can be reproduced by running `Examples.R`. This file also serves as a practical guide for using the functions described above.

```r
source("Functions.R")
source("Examples.R")
```

## BioModels sweep

The benchmark added in the revised manuscript searches for suitable preorders on 205 reaction networks obtained from BioModels. The inputs, recorded results and instructions for repeating the computations are in [Sweep/README.md](Sweep/README.md).
