#!/bin/bash
# Usage: stochasticbarrier pwc

if [ "$1" != "pwc" ]; then
    echo "Usage: $0 pwc"
    exit 1
fi

# -------------------------------
# Julia SOS and PWC experiments
# -------------------------------
JULIA_PROJECT="/StochasticBarrierFunctions"
echo "Starting Julia SOS vs PWC experiments ..."

julia --project="$JULIA_PROJECT" <<'EOF'
# Define the YAML files for SOS
yaml_files = [

    # ---------------------SOS Comparison Systems---------------------#
    "benchmarks/linear/systems/contraction2/SOS/sos_deg4.yaml",
    "benchmarks/linear/systems/contraction/SOS/sos_deg8.yaml"
    # "benchmarks/linear/systems/contraction/SOS/sos_deg20.yaml",
    # "benchmarks/linear/systems/contraction/SOS/sos_deg30.yaml"

]

# Include the barrier synthesis once
include("benchmarks/barrier_synthesis.jl")

# Run all Julia experiments
for yaml_file in yaml_files
    println("Running Julia experiment: ", yaml_file)
    barrier_synthesis(yaml_file)
end

println("All Julia PWC vs SOS experiments completed!")
EOF
