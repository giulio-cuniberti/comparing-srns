# BioModels sweep

This folder contains the inputs, code and results for the computation-time table in the paper. The search checks the sufficient conditions of Theorem 4.1 on 205 reaction networks obtained from [BioModels](https://www.ebi.ac.uk/biomodels/), with at most ten species. It finds at least one non-trivial comparison for 133 networks; an unsuccessful search does not rule out other comparisons.

## Files

- `Run.R`: runs the searches or summarizes saved results, using `../Functions.R`.
- `Networks.R`: source and product matrices, BioModels identifiers and model names.
- `Results.csv`: the results reported in the paper.

The source and product matrices `S` and `P` have one row per reaction and one column per species, as described in the [main README](../README.md). The sweep uses these matrices; it does not include or check the original kinetic laws. Applying the comparisons requires the kinetic assumptions of Theorem 4.1.

## Running the computations

Install R and the packages `rcdd` and `gmp` as described in the main README. Run these commands in a terminal **from the repository folder**.

Display the recorded table without repeating the searches:

```sh
Rscript --vanilla Sweep/Run.R --summary
```

Try the networks with at most four species:

```sh
Rscript --vanilla Sweep/Run.R --max-species 4 --output Sweep/Results-test.csv
```

Repeat the full sweep:

```sh
Rscript --vanilla Sweep/Run.R
```

The default is one worker for networks with at most five species and seven for larger networks. Windows uses one worker. Add `--workers 1` to run serially on any system.

New results are saved in `Sweep/Results-new.csv`; `Results.csv` is preserved. Repeating a command resumes the sweep, skipping successfully completed networks. Keep the inputs, code and worker settings unchanged when resuming; choose a new `--output` filename to start again. Relative paths supplied with `--output` start from the working directory. To summarize new results, run:

```sh
Rscript --vanilla Sweep/Run.R --summary Sweep/Results-new.csv
```

## Search and timings

For each species, the search considers either direction of comparison, equality, or no comparison. Choices differing only by reversing all inequalities are counted once, giving `(4^d + 2^d)/2` candidate matrices for `d` species. All candidates are checked. The two trivial comparisons (all species equal, or no species compared) are excluded from the reported successes.

The benchmark checked **3,243,525 candidates** on 25 September 2026 in approximately **1 hour 49 minutes**, using the default worker settings on an **Apple M3 MacBook Air with 16 GB RAM**. The software was macOS 27.0, R 4.4.1, `rcdd` 1.6-1, `gmp` 0.7-5.1 and GMP 6.3.0.

Per-network times are rounded to milliseconds and cover candidate checks, excluding preliminary setup and final processing. Timings depend on the machine and workload.

## Columns in `Results.csv`

Comparisons obtained by reversing all inequalities are counted once.

| Column | Meaning |
|---|---|
| `source`, `id`, `name` | Input group, BioModels identifier and model name. |
| `d`, `n` | Number of species and reactions. |
| `candidates` | Number of candidate matrices checked. |
| `lp` | Number of cone computations used to check the conditions on `A` and `B`. |
| `wall_s` | Elapsed time for candidate checks, in seconds. |
| `cpu_s` | CPU time of the main process and its workers, in seconds. |
| `cores` | Number of workers used. |
| `orderable` | `1` if a non-trivial preorder or equivalence was found; `0` otherwise. |
| `n_preorders` | Number of distinct successful comparisons containing a `+` or `-` species inequality. |
| `n_equivalences` | Number of distinct successful comparisons using only equality and unconstrained species. |
| `raw_hits` | Number of distinct successful comparisons before removing the two trivial ones. |
| `status` | `ok` means that the search completed successfully. |
