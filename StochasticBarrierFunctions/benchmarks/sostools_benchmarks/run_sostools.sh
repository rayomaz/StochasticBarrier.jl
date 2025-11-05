#!/bin/bash
# =================================================
# run_sostools.sh - Run all SOSTOOLS MATLAB benchmarks
# =================================================

# Exit immediately if any command fails
set -e

# Base directory inside container
BASE_DIR="/app"

# List of benchmarks to run
benchmarks=(
    "linear/contraction.m"
)

# Run each benchmark in MATLAB
for bench in "${benchmarks[@]}"; do
    echo "------------------------------------------------------------"
    echo "Running SOSTOOLS benchmark: $bench"
    echo "------------------------------------------------------------"

    matlab -batch "run('${BASE_DIR}/${bench}')"

    echo "------------------------------------------------------------"
    echo "Completed benchmark: $bench"
done
