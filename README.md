# Comparing SRNs

This repository contains two R scripts:
- `Functions.R`
- `Examples.R`

## Installation and Usage

After sourcing `Functions.R`, you can use its two main functions:

### `checkOrdering(source.complexes, product.complexes, preorder.matrix)`
Checks whether Theorem 4.1 can be applied to the reaction network defined by the first two arguments, given the preorder specified by the third argument.

### `findOrderings(source.complexes, product.complexes)`
Finds all matrix preorders that enable the application of Theorem 4.1 to the reaction network defined by the arguments.

## Input Format

- `source.complexes` and `product.complexes` are `n × d` matrices, where:
  - `n` = number of reactions
  - `d` = number of species
- Row `r` in `source.complexes` contains the stoichiometric coefficients of the **source complex** for reaction `r`
- Row `r` in `product.complexes` contains the stoichiometric coefficients of the **product complex** for the same reaction `r`
- `preorder.matrix` is the matrix *M* from Theorem 4.1

## Output Description

[... DESCRIZIONE DEGLI OUTPUT QUI ...]

## Examples

All examples presented in the paper can be reproduced using `Examples.R`. This file also serves as a practical guide for using the functions described above.
