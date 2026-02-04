# Comparing SRNs

This repository contains the R implementation for the paper **"Stochastic ordering tools for continuous-time Markov chains and applications to reaction network models"** by Daniele Cappelletti, Giulio Cuniberti and Paola Siri. The code consists of two main scripts:
- `Examples.R`
- `Functions.R`

## Installation and usage

After sourcing `Functions.R`, you can use its two main functions:

#### `checkOrdering(source.complexes, product.complexes, preorder.matrix)`
checks whether Theorem 4.1 can be applied to the reaction network defined by the first two arguments, given the preorder specified by the third argument.

#### `findOrderings(source.complexes, product.complexes)`
finds all rate inequalities and associated preorders on species counts that enable the application of Theorem 4.1 to the reaction network defined by the arguments.

## Input format

- `source.complexes` and `product.complexes` are `n × d` matrices, where:
  - `n` = number of reactions
  - `d` = number of species
  - row `r` in `source.complexes` represents the **source complex** of reaction `r`
  - row `r` in `product.complexes` represents the **product complex** of the same reaction `r`
  - column `s` refers to the same species `s`, in both matrices

- `preorder.matrix` is the matrix *M* from Theorem 4.1

## Output description

**`checkOrdering(source.complexes, product.complexes, preorder.matrix)`** returns a list of 4 objects.
- `result`: 1 or 0, depending on whether the hypotheses of Theorem 4.1 are satisfied or not
- `rate.inequalities`: vector whose component `r` corresponds to reaction `r` and can take the values
  - `"="` if the rate constants are equal in *X* and *Y*
  - `"+"` if the rate constant of *Y* is greater or equal than that of *X*
  - `"-"` if the rate constant of *Y* is less or equal than that of *X*
  - `"?"` if the rate constants of *X* and *Y* are not compared
- `species.inequalities`: vector whose component `s` corresponds to species `s` and can take the values
  - `"="` if the molecular counts are equal in *X* and *Y*
  - `"+"` if the molecular count of *Y* is greater or equal than that of *X*
  - `"-"` if the molecular count of *Y* is less or equal than that of *X*
  - `"?"` if the molecular counts of *X* and *Y* are not compared
- `first.fail`: when `result = 0`, it specifies the first hypothesis of Theorem 4.1 that failed to be satisfied

**`findOrderings(source.complexes, product.complexes)`** returns a list of 2 objects.
- `preorders`: list of all pairs of `rate.inequalities` and `species.inequalities` (as described for the previous function) that define an admissible preordering structure
  and for which at least one component of `species.inequalities` is either `"+"` or `"-"`
- `equivalences`: list of all pairs of `rate.inequalities` and `species.inequalities` (as described for the previous function) that define the remaining admissible preordering (equivalence) structures

## Examples

All examples presented in the paper can be reproduced by running `Examples.R`. This file also serves as a practical guide for using the functions described above.
