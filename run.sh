#!/usr/bin/env bash
# run.sh -- launch (or resume) the Arm A full run in the background.
# Usage: ./run.sh [cores]   (default: all cores minus 2)
# Resumes from armA-rct/results/checkpoint_200reps.rds if present.
set -euo pipefail
cd "$(dirname "$0")/armA-rct"
CORES="${1:-$(( $(nproc) - 2 ))}"
# pin BLAS/OpenMP threads BEFORE R starts: without this the forked workers
# oversubscribe and the run takes about 30 times longer
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p results
nohup Rscript sim-run.R 200 "$CORES" > results/sim-run.log 2>&1 &
echo "PID $! on $CORES cores. Follow with: tail -c 300 armA-rct/results/sim-run.log"
