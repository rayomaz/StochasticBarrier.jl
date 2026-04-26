#!/bin/bash
set -e

echo "------------------------------------------------------------"
echo "Setting up SOSTOOLS with MOSEK"
echo "------------------------------------------------------------"

# === Clone SOSTOOLS ===
if [ ! -d "SOSTOOLS" ]; then
    git clone https://github.com/oxfordcontrol/SOSTOOLS.git
else
    echo "SOSTOOLS already exists — skipping clone."
fi

# === Download MOSEK ===
if [ ! -d "mosek" ]; then
    echo "------------------------------------------------------------"
    echo "Downloading MOSEK"
    echo "------------------------------------------------------------"

    wget https://download.mosek.com/stable/10.1.21/mosektoolslinux64x86.tar.bz2
    tar -xjf mosektoolslinux64x86.tar.bz2
else
    echo "MOSEK already exists — skipping download."
fi

# === Copy MOSEK license ===
echo "------------------------------------------------------------"
echo "Copying MOSEK license"
echo "------------------------------------------------------------"

mkdir -p mosek
cp -f ../mosek/mosek.lic mosek/mosek.lic

# === MATLAB setup ===
echo "------------------------------------------------------------"
echo "Configuring MATLAB"
echo "------------------------------------------------------------"

matlab -nodisplay -nosplash -nodesktop <<'EOF'
addpath('SOSTOOLS');
addsostools;

addpath(genpath('mosek'));

disp('Checking MOSEK...');
mosekopt('version')

savepath;
exit
EOF

# === Run benchmarks ===
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

    matlab -nodisplay -nosplash -nodesktop <<EOF
addpath('SOSTOOLS');
addpath(genpath('mosek'));

solver_opt = sosoptions;
solver_opt.solver = 'mosek';

run('$EXP');
exit
EOF

done

echo ""
echo "------------------------------------------------------------"
echo "All SOSTOOLS benchmarks completed with MOSEK."
echo "------------------------------------------------------------"