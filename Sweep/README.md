# BioModels sweep

This folder contains the computations accompanying the runtime table in the revised manuscript. The search tests the sufficient conditions of Theorem 4.1 on 205 reaction networks with at most ten species, finding at least one non-trivial preorder or equivalence for 133 of them. An unsuccessful search does not rule out other stochastic comparisons.

## Files

- `Run.R`: runs the searches using the functions in `../Functions.R`, or summarizes saved results.
- `Networks.R`: source and product matrices for the 205 networks, with their BioModels identifiers and names.
- `Results.csv`: the recorded results underlying the manuscript's table.

Each input has one row per reaction and one column per species. Its matrices `S` and `P` have the same meaning as `source.complexes` and `product.complexes` in the [main README](../README.md). The search uses these matrices, not numerical rate constants or propensity formulas.

The inputs are the matrices obtained from [BioModels](https://www.ebi.ac.uk/biomodels/) and used in the benchmark, preserving their row and column order; each identifier occurs once. They are provided so that the same matrix computations can be repeated without downloading or re-extracting the database. These are reaction-network data: the original kinetic laws are not included or checked by this script. The labels `HARV`, `MM`, `HILL` and `BM` retain the source groups of the extraction.

## Inspecting one network

From the repository folder, the first network can be examined in R without running the full sweep:

```r
source("Functions.R")
source("Sweep/Networks.R")
z <- networks[[1]]
findOrderings(z$S, z$P)
```

## Running the computations

Install R and the two packages described in the main README. Run the following commands in a terminal from the repository folder.

To display the manuscript's summary table from the recorded results, without repeating the searches:

```sh
Rscript --vanilla Sweep/Run.R --summary
```

For a small trial using only networks with at most four species:

```sh
Rscript --vanilla Sweep/Run.R --workers 1 --max-species 4 --output Sweep/Results-test.csv
```

To repeat the full sweep:

```sh
Rscript --vanilla Sweep/Run.R
```

The full sweep uses one worker for networks with at most five species and seven workers for larger networks. On Windows it runs serially. Use `--workers 1` to run serially on any system.

New results are saved as `Results-new.csv` in this folder, leaving `Results.csv` unchanged. Repeating a command resumes its output file, skipping networks already completed successfully; the worker settings must agree. Choose a new `--output` filename to repeat all searches. A supplied output path is relative to the terminal's working directory. To summarize a new file, use `Rscript --vanilla Sweep/Run.R --summary Sweep/Results-new.csv`.

## Search and timings

For each species, the search considers either direction of comparison, equality, or no comparison. Choices differing only by reversing all inequalities are counted once. Thus a network with `d` species has `(4^d + 2^d)/2` candidate matrices. Every candidate is checked, even after a successful comparison has been found. The two trivial comparisons are included among the candidates but excluded from the reported successes.

The recorded run checked **3,243,525 candidates** on 25 September 2026. The full run took approximately **1 hour 49 minutes** (6,521 seconds measured by the driver), on a MacBook Air with an **Apple M3 processor and 16 GB RAM**. It used macOS 27.0, R 4.4.1, `rcdd` 1.6-1, `gmp` 0.7-5.1 and GMP 6.3.0, with the worker policy described above.

The per-network times in `Results.csv` measure the checks of all candidate matrices, excluding preliminary setup and final merging of results. Their sum is 6,449.490 seconds; it is therefore smaller than the duration of the full run. Times are recorded to the nearest millisecond, so a displayed zero represents a very short computation. Results should be reproducible; timings depend on the machine and workload.

## Columns in `Results.csv`

| Column | Meaning |
|---|---|
| `source`, `id`, `name` | Extraction group, BioModels identifier and model name. |
| `d`, `n` | Number of species and reactions in the supplied matrices. |
| `candidates` | Number of candidate matrices checked. |
| `lp` | Number of calls computing a cone description for the conditions on `A` and `B`; the column retains its original name and does not count individual linear programs. |
| `wall_s` | Elapsed search time, in seconds, with the scope described above. |
| `cpu_s` | Recorded CPU time of the main process and its workers, in seconds. |
| `cores` | Number of workers used. |
| `orderable` | `1` if at least one non-trivial preorder or equivalence was found, otherwise `0`. |
| `n_preorders` | Number of distinct successful comparisons containing a `+` or `-` species inequality. |
| `n_equivalences` | Number of distinct successful comparisons using only equality and unconstrained species. |
| `raw_hits` | Number of distinct successful comparisons before removing the two trivial ones. |
| `status` | `ok` confirms that the search completed successfully. |

The two comparison counts use the same convention as `findOrderings`: reversed comparisons are identified. A successful search is one for which `n_preorders > 0` or `n_equivalences > 0`.
