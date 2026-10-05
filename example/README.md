# scLongTree example

This directory contains a synthetic two-timepoint example for testing
scLongTree.

## Input data

The example is the `default_1/rep1` synthetic dataset used in the scLongTree
simulation experiments.

- `input.D.csv` - combined D matrix: 398 cells x 74 mutations
- `t1.D.csv` - timepoint 1 matrix: 300 cells x 74 mutations
- `t2.D.csv` - timepoint 2 matrix: 98 cells x 74 mutations
- `cell_timepoints.csv` - assignment of rows in `input.D.csv` to t1 and t2

`t1.D.csv` followed by `t2.D.csv` exactly reconstructs `input.D.csv`.

All data in this example are simulated.

## Precomputed BnpC results

`bnpc_runs/` contains five precomputed BnpC runs for each timepoint:

    bnpc_runs/
      m1/t1/
      m1/t2/
      ...
      m5/t1/
      m5/t2/

These results are included so users can test scLongTree without rerunning
BnpC.

Each BnpC run contains:

- `assignment.txt`
- `errors.txt`
- `genotypes_posterior_mean.tsv`
- `genotypes_cont_posterior_mean.tsv`

## Quick test

From the repository root:

    bash run.sh \
      --python /path/to/python \
      --bnpc-results example/bnpc_runs \
      --D example/input.D.csv \
      --cells example/cell_timepoints.csv \
      --out example_test \
      --sample default_1_rep1 \
      --k 0 \
      --p 0

The resulting tree and plot are written under `TreeCSVs/` and `plots/`.

See the main repository `README.md` for the full workflow, including how to
rerun BnpC from the timepoint matrices.
