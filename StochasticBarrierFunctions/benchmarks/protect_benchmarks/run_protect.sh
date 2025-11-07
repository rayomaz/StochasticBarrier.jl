#!/bin/bash
# ============================================================
# run_protect.sh - Benchmark PRoTECT experiments
# ============================================================

set -e

# --- Ensure PYTHONPATH includes PRoTECT src/ ---
export PYTHONPATH="/app/PRoTECT/src:$PYTHONPATH"

# Change to PRoTECT root
cd /app/PRoTECT

# --- Run linear experiments as modules ---
python3 -m experiments.linear.contractive
python3 -m experiments.linear.twotank
python3 -m experiments.linear.room
python3 -m experiments.linear.quadrotor

# --- Run polynomial experiments as modules ---
python3 -m experiments.polynomial.thermostat
python3 -m experiments.polynomial.oscillator
