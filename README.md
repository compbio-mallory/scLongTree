# scLongTree

scLongTree is a computational method for inferring longitudinal clonal evolution trees from single-cell DNA sequencing (scDNA-seq) data collected at multiple timepoints.

scLongTree uses clustering and genotype estimates from multiple BnpC runs and integrates information across timepoints to infer a longitudinal clonal tree. The method can additionally model back mutations (mutation losses) and recurrent/parallel mutation events.

## Repository contents

This repository contains the scLongTree source code and a synthetic example dataset for testing the software.

```text
scLongTree/
├── algorithm/                  # scLongTree source code
├── example/
│   ├── input.D.csv             # combined cell-by-mutation matrix
│   ├── t1.D.csv                # timepoint 1 matrix
│   ├── t2.D.csv                # timepoint 2 matrix
│   ├── cell_timepoints.csv     # cell-to-timepoint assignments
│   └── bnpc_runs/              # precomputed BnpC results
├── run.sh                      # wrapper for running the pipeline
└── README.md
```

## Simulator
Simulation utilities used to generate synthetic datasets are available on the
[`simulator` branch](https://github.com/compbio-mallory/scLongTree/tree/simulator).


## Requirements

scLongTree requires Python 3. The current example has been tested with Python 3.8.

BnpC is used to obtain cell clusters and clone genotypes at individual timepoints.

Precomputed BnpC outputs are included with the example, so BnpC does **not** need to be rerun simply to test scLongTree.

Users who wish to perform the complete workflow from raw timepoint matrices will additionally need a working BnpC installation.

## Example dataset

The `example/` directory contains a synthetic two-timepoint dataset from the scLongTree simulation experiments.

Five precomputed BnpC runs are supplied for each timepoint:

```text
example/bnpc_runs/
├── m1/
│   ├── t1/
│   └── t2/
├── m2/
├── m3/
├── m4/
└── m5/
```

Each run contains the BnpC assignments, error-rate estimates, and posterior genotype estimates required by scLongTree.

## Quick start

The fastest way to test scLongTree is to use the bundled BnpC results.

From the repository root:

```bash
bash run.sh \
  --python /path/to/python \
  --bnpc-results example/bnpc_runs \
  --D example/input.D.csv \
  --cells example/cell_timepoints.csv \
  --out example_test \
  --sample default_1_rep1 \
  --k 0 \
  --p 0
```

## Running BnpC from scratch

The bundled BnpC outputs are provided only to make testing faster. Users can instead run the complete workflow by supplying BnpC's `run_BnpC.py` script.

```bash
bash run.sh \
  --python /path/to/python \
  --bnpc /path/to/run_BnpC.py \
  --t1 example/t1.D.csv \
  --t2 example/t2.D.csv \
  --D example/input.D.csv \
  --cells example/cell_timepoints.csv \
  --out example_full \
  --sample default_1_rep1 \
  --k 0 \
  --p 0
```

By default, the wrapper performs five BnpC runs at each timepoint.

The number of BnpC runs and the number executed concurrently can be changed using:

```text
--bnpc-runs N
--parallel-bnpc N
```

## scLongTree parameters

The main evolutionary-event parameters are:

### Back mutations

```text
--k N
```

sets the maximum number of allowed mutation-loss/back-mutation events.

For example:

```text
--k 0
```

disallows back mutations.

### Parallel mutations

```text
--p N
```

sets the maximum number of additional parallel occurrences allowed for a mutation.

The total number of allowed occurrences for a mutation is:

```text
1 + p
```

Thus:

```text
--p 0
```

allows each mutation to occur once, while:

```text
--p 1
```

allows one additional parallel occurrence.

## Running scLongTree directly

scLongTree can also be called directly using `algorithm/selectBestTree.py`.

For example:

```bash
python algorithm/selectBestTree.py \
  -m 5 \
  -t t1 t2 \
  -loc example/bnpc_runs \
  -cells example/cell_timepoints.csv \
  -D example/input.D.csv \
  -k 0 \
  -p 0 \
  -op example_test \
  -sample default_1_rep1
```


- `-m`: number of BnpC runs
- `-t`: ordered timepoint labels
- `-loc`: directory containing the BnpC results
- `-cells`: file assigning cells to timepoints
- `-D`: combined cell-by-mutation matrix
- `-k`: maximum number of back mutations/losses
- `-p`: maximum number of additional parallel mutation occurrences
- `-op`: output prefix
- `-sample`: sample name used when generating the tree visualization

The underlying scLongTree program supports multiple timepoints through the `-t` argument. 

## Input format

### D matrix

The D matrix is a tab-delimited cell-by-mutation matrix.

Rows represent cells and columns represent mutations.

Example:

```text
0    1    0    3
0    1    1    0
1    1    0    0
```

### Cell-timepoint file

The cell-timepoint file contains one line per timepoint. Cell indices correspond to rows of the combined D matrix.

For example:

```text
t1    0;1;2;3
t2    4;5;6
```
