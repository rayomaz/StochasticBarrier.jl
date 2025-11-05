#!/bin/bash
set -e

echo "------------------------------------------------------------"
echo "Setting up SOSTOOLS and SDPT3 environment"
echo "------------------------------------------------------------"

# === Clone SOSTOOLS ===
if [ ! -d "SOSTOOLS" ]; then
    git clone https://github.com/oxfordcontrol/SOSTOOLS.git
else
    echo "SOSTOOLS already exists — skipping clone."
fi

# === Clone SDPT3 ===
if [ ! -d "SDPT3" ]; then
    git clone https://github.com/sqlp/sdpt3.git SDPT3
else
    echo "SDPT3 already exists — skipping clone."
fi

# === Run MATLAB setup ===
echo "------------------------------------------------------------"
echo "Installing SOSTOOLS and SDPT3"
echo "------------------------------------------------------------"
matlab -batch "addpath('SOSTOOLS'); addsostools; addpath('SDPT3'); install_sdpt3; savepath; exit"

# === Run all benchmark experiments ===
echo "------------------------------------------------------------"
echo "Running SOSTOOLS Benchmarks"
echo "------------------------------------------------------------"

EXPERIMENTS=(
    "linear/contraction.m"
    "linear/twotank.m"
    "linear/room.m"
)

for EXP in "${EXPERIMENTS[@]}"; do
    echo ""
    echo "------------------------------------------------------------"
    echo "Running experiment: $EXP"
    echo "------------------------------------------------------------"
    matlab -batch "run('$EXP'); exit"
done

echo ""
echo "------------------------------------------------------------"
echo "All SOSTOOLS benchmarks completed."
echo "------------------------------------------------------------"
