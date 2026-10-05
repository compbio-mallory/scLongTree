#!/usr/bin/env bash

set -Eeuo pipefail

# ============================================================
# scLongTree pipeline
#
# Two modes are supported:
#
# 1. QUICK MODE
#    Use previously generated BnpC results.
#    This is the recommended way to test the bundled example.
#
# 2. FULL MODE
#    Run BnpC first, then run scLongTree.
#
# The current wrapper handles two observed timepoints: t1 and t2.
# scLongTree itself supports multiple timepoints through selectBestTree.py.
# ============================================================


usage() {
cat <<'USAGE'

Usage
=====

QUICK MODE
----------
Use precomputed BnpC results:

  bash run.sh \
    --python /path/to/python \
    --bnpc-results example/bnpc_runs \
    --D example/input.D.csv \
    --cells example/cell_timepoints.csv \
    --out example_test \
    --sample example

FULL MODE
---------
Run BnpC first and then scLongTree:

  bash run.sh \
    --python /path/to/python \
    --bnpc /path/to/run_BnpC.py \
    --t1 example/t1.D.csv \
    --t2 example/t2.D.csv \
    --D example/input.D.csv \
    --cells example/cell_timepoints.csv \
    --out example_full \
    --sample example

Common required arguments
-------------------------

  --python PATH
      Python executable used for BnpC/scLongTree.

  --D FILE
      Combined cell-by-mutation D matrix.

  --cells FILE
      Cell-to-timepoint assignment file used by scLongTree.

  --out PREFIX
      Output prefix relative to the scLongTree repository.

      For example:

          --out example_test

      produces outputs under:

          output/example_test/
          TreeCSVs/example_test*
          plots/example_test*

Quick-mode argument
-------------------

  --bnpc-results DIR
      Directory containing existing BnpC runs in the form:

          DIR/
            m1/t1/
            m1/t2/
            ...
            m5/t1/
            m5/t2/

Full-mode arguments
-------------------

  --bnpc FILE
      BnpC run_BnpC.py script.

  --t1 FILE
      D matrix for timepoint t1.

  --t2 FILE
      D matrix for timepoint t2.

Optional arguments
------------------

  --sample NAME
      Sample name used for plotting.
      Default: sample

  --bnpc-runs N
      Number of BnpC runs.
      Default: 5

  --parallel-bnpc N
      Maximum number of BnpC runs executed simultaneously.
      Default: 5

  --k N
      Maximum number of back mutations/losses.
      Default: 0

  --p N
      Maximum number of extra parallel mutation occurrences.
      Total allowed occurrences for a locus = 1 + P.
      Default: 0

  -h, --help
      Show this help message.

USAGE
}


# ============================================================
# Defaults
# ============================================================

PYTHON=""

BNPC=""
BNPC_RESULTS=""

T1=""
T2=""

FULL_D=""
CELLS=""

OUT_PREFIX=""
SAMPLE="sample"

BNPC_RUNS=5
PARALLEL_BNPC=5

K=0
P=0


# ============================================================
# Parse arguments
# ============================================================

while [[ $# -gt 0 ]]; do

    case "$1" in

        --python)
            PYTHON="$2"
            shift 2
            ;;

        --bnpc)
            BNPC="$2"
            shift 2
            ;;

        --bnpc-results)
            BNPC_RESULTS="$2"
            shift 2
            ;;

        --t1)
            T1="$2"
            shift 2
            ;;

        --t2)
            T2="$2"
            shift 2
            ;;

        --D)
            FULL_D="$2"
            shift 2
            ;;

        --cells)
            CELLS="$2"
            shift 2
            ;;

        --out)
            OUT_PREFIX="$2"
            shift 2
            ;;

        --sample)
            SAMPLE="$2"
            shift 2
            ;;

        --bnpc-runs)
            BNPC_RUNS="$2"
            shift 2
            ;;

        --parallel-bnpc)
            PARALLEL_BNPC="$2"
            shift 2
            ;;

        --k)
            K="$2"
            shift 2
            ;;

        --p)
            P="$2"
            shift 2
            ;;

        -h|--help)
            usage
            exit 0
            ;;

        *)
            echo "ERROR: unknown argument: $1"
            echo
            usage
            exit 1
            ;;

    esac

done


# ============================================================
# Common validation
# ============================================================

for var in PYTHON FULL_D CELLS OUT_PREFIX; do

    if [[ -z "${!var}" ]]; then

        echo "ERROR: required argument is missing: $var"
        echo
        usage
        exit 1

    fi

done


# Python can either be an absolute executable or a command in PATH.

if [[ ! -x "$PYTHON" ]] \
   && ! command -v "$PYTHON" >/dev/null 2>&1; then

    echo "ERROR: Python executable not found:"
    echo "  $PYTHON"
    exit 1

fi


for f in "$FULL_D" "$CELLS"; do

    if [[ ! -s "$f" ]]; then

        echo "ERROR: required input file is missing or empty:"
        echo "  $f"
        exit 1

    fi

done


if (( BNPC_RUNS < 1 )); then

    echo "ERROR: --bnpc-runs must be >= 1"
    exit 1

fi


if (( PARALLEL_BNPC < 1 )); then

    echo "ERROR: --parallel-bnpc must be >= 1"
    exit 1

fi


if (( PARALLEL_BNPC > BNPC_RUNS )); then
    PARALLEL_BNPC="$BNPC_RUNS"
fi


# Keep scLongTree's output prefix repository-relative.

if [[ "$OUT_PREFIX" = /* ]]; then

    echo "ERROR: --out must be a relative output prefix."
    echo
    echo "Example:"
    echo "  --out example_test"
    exit 1

fi


# ============================================================
# Repository/output locations
# ============================================================

SCRIPT_DIR="$(
    cd "$(dirname "${BASH_SOURCE[0]}")"
    pwd
)"

RUN_DIR="$SCRIPT_DIR/output/$OUT_PREFIX"
LOG_DIR="$RUN_DIR/logs"

mkdir -p "$RUN_DIR" "$LOG_DIR"

mkdir -p \
    "$SCRIPT_DIR/TreeCSVs/$(dirname "$OUT_PREFIX")" \
    "$SCRIPT_DIR/plots/$(dirname "$OUT_PREFIX")"


# ============================================================
# Select mode
# ============================================================

if [[ -n "$BNPC_RESULTS" ]]; then

    MODE="precomputed"

    if [[ -n "$BNPC" || -n "$T1" || -n "$T2" ]]; then

        echo "ERROR:"
        echo "  --bnpc-results cannot be combined with"
        echo "  --bnpc, --t1, or --t2."
        exit 1

    fi

    if [[ ! -d "$BNPC_RESULTS" ]]; then

        echo "ERROR: BnpC results directory does not exist:"
        echo "  $BNPC_RESULTS"
        exit 1

    fi

    BNPC_LOC="$(
        cd "$BNPC_RESULTS"
        pwd
    )"

else

    MODE="run_bnpc"

    for var in BNPC T1 T2; do

        if [[ -z "${!var}" ]]; then

            echo "ERROR:"
            echo "  Full mode requires --bnpc, --t1, and --t2."
            exit 1

        fi

    done

    for f in "$BNPC" "$T1" "$T2"; do

        if [[ ! -s "$f" ]]; then

            echo "ERROR: required full-mode file is missing or empty:"
            echo "  $f"
            exit 1

        fi

    done

    BNPC_LOC="$RUN_DIR/bnpc"

    mkdir -p "$BNPC_LOC"

fi


# ============================================================
# Save run configuration
# ============================================================

CONFIG="$RUN_DIR/run_config.txt"

cat > "$CONFIG" <<CONFIG
mode=$MODE
python=$PYTHON
D=$FULL_D
cells=$CELLS
out_prefix=$OUT_PREFIX
sample=$SAMPLE
bnpc_runs=$BNPC_RUNS
parallel_bnpc=$PARALLEL_BNPC
k=$K
p=$P
bnpc_location=$BNPC_LOC
date=$(date)
host=$(hostname)
CONFIG


echo
echo "============================================================"
echo "scLongTree"
echo "============================================================"
echo "Mode:             $MODE"
echo "Sample:           $SAMPLE"
echo "BnpC runs:        $BNPC_RUNS"
echo "k:                $K"
echo "p:                $P"
echo "BnpC location:    $BNPC_LOC"
echo "Output prefix:    $OUT_PREFIX"
echo "============================================================"
echo


# ============================================================
# Full mode: run BnpC
# ============================================================

run_one_bnpc() {

    local run="$1"
    local tp="$2"
    local input="$3"

    local odir="$BNPC_LOC/m${run}/${tp}"

    mkdir -p "$odir"

    echo "[BnpC] starting m${run}/${tp}"

    set +e

    "$PYTHON" "$BNPC" \
        "$input" \
        -t \
        -o "$odir" \
        -np \
        > "$odir/bnpc.stdout.log" \
        2> "$odir/bnpc.stderr.log"

    local status=$?

    set -e

    echo "$status" > "$odir/exit_status.txt"

    if [[ "$status" -ne 0 ]]; then

        echo "ERROR: BnpC failed for m${run}/${tp}"
        echo "See:"
        echo "  $odir/bnpc.stderr.log"

        return "$status"

    fi

    echo "[BnpC] finished m${run}/${tp}"

}


run_bnpc_timepoint() {

    local tp="$1"
    local input="$2"

    echo
    echo "Running BnpC for $tp"
    echo "Input: $input"
    echo

    local pids=()
    local labels=()
    local failed=0

    for run in $(seq 1 "$BNPC_RUNS"); do

        run_one_bnpc \
            "$run" \
            "$tp" \
            "$input" &

        pids+=("$!")
        labels+=("m${run}/${tp}")


        if [[ "${#pids[@]}" -ge "$PARALLEL_BNPC" ]]; then

            for i in "${!pids[@]}"; do

                if ! wait "${pids[$i]}"; then

                    echo "ERROR: ${labels[$i]} failed."
                    failed=1

                fi

            done

            pids=()
            labels=()

        fi

    done


    # Wait for final partial batch.

    for i in "${!pids[@]}"; do

        if ! wait "${pids[$i]}"; then

            echo "ERROR: ${labels[$i]} failed."
            failed=1

        fi

    done


    if [[ "$failed" -ne 0 ]]; then

        echo "ERROR: one or more BnpC runs failed for $tp."
        exit 1

    fi

}


if [[ "$MODE" == "run_bnpc" ]]; then

    echo "============================================================"
    echo "Running BnpC"
    echo "============================================================"

    # Run the two timepoints sequentially.
    # Multiple BnpC replicates within each timepoint may run in parallel.

    run_bnpc_timepoint "t1" "$T1"
    run_bnpc_timepoint "t2" "$T2"

fi


# ============================================================
# Verify BnpC results
# ============================================================

echo
echo "============================================================"
echo "Verifying BnpC results"
echo "============================================================"

for run in $(seq 1 "$BNPC_RUNS"); do

    for tp in t1 t2; do

        d="$BNPC_LOC/m${run}/${tp}"

        for required in \
            assignment.txt \
            errors.txt \
            genotypes_posterior_mean.tsv \
            genotypes_cont_posterior_mean.tsv
        do

            if [[ ! -s "$d/$required" ]]; then

                echo "ERROR: missing BnpC result:"
                echo "  $d/$required"
                exit 1

            fi

        done

        echo "PASS: m${run}/${tp}"

    done

done


# ============================================================
# Run scLongTree
# ============================================================

echo
echo "============================================================"
echo "Running scLongTree"
echo "============================================================"

SCLT_STDOUT="$LOG_DIR/sclongtree.stdout.log"
SCLT_STDERR="$LOG_DIR/sclongtree.stderr.log"

cd "$SCRIPT_DIR"

set +e

"$PYTHON" \
    algorithm/selectBestTree.py \
    -m "$BNPC_RUNS" \
    -t t1 t2 \
    -loc "$BNPC_LOC" \
    -cells "$CELLS" \
    -D "$FULL_D" \
    -k "$K" \
    -p "$P" \
    -op "$OUT_PREFIX" \
    -sample "$SAMPLE" \
    > "$SCLT_STDOUT" \
    2> "$SCLT_STDERR"

SCLT_STATUS=$?

set -e


if [[ "$SCLT_STATUS" -ne 0 ]]; then

    echo
    echo "ERROR: scLongTree failed with exit status $SCLT_STATUS"
    echo
    echo "stderr:"
    echo "  $SCLT_STDERR"
    echo
    echo "Last 30 stderr lines:"
    tail -30 "$SCLT_STDERR" || true

    exit "$SCLT_STATUS"

fi


# ============================================================
# Success
# ============================================================

echo
echo "============================================================"
echo "scLongTree completed successfully"
echo "============================================================"
echo
echo "Run information:"
echo "  $RUN_DIR"
echo
echo "scLongTree stdout:"
echo "  $SCLT_STDOUT"
echo
echo "Tree output prefix:"
echo "  $SCRIPT_DIR/TreeCSVs/$OUT_PREFIX"
echo
echo "Plot output prefix:"
echo "  $SCRIPT_DIR/plots/$OUT_PREFIX"
echo
echo "============================================================"
