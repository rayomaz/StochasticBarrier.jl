#!/bin/bash
# ============================================================
# run_protect.sh - Benchmark PRoTECT experiments
# ============================================================

set -e

# --- Ensure PYTHONPATH includes PRoTECT src/ ---
export PYTHONPATH="/app/PRoTECT/src:$PYTHONPATH"

# --- Run linear and polynomial experiments ---
python3 /app/PRoTECT/experiments/linear/contractive.py
python3 /app/PRoTECT/experiments/linear/twotank.py
python3 /app/PRoTECT/experiments/linear/room.py
python3 /app/PRoTECT/experiments/linear/quadrotor.py

python3 /app/PRoTECT/experiments/polynomial/thermostat.py
python3 /app/PRoTECT/experiments/polynomial/oscillator.py
