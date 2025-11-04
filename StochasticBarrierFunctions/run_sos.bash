#!/bin/bash
# Usage: stochasticbarrier sos
# Single Julia call runs all SOS experiments

if [ "$1" != "sos" ]; then
    echo "Usage: $0 sos"
    exit 1
fi

JULIA_PROJECT="/StochasticBarrierFunctions"
echo "Starting SOS experiments ..."

julia --project="$JULIA_PROJECT" <<'EOF'
# Define the YAML files
yaml_files = [
    "benchmarks/linear/systems/contraction/SOS/sos_deg2.yaml",
    "benchmarks/linear/systems/contraction/SOS/sos_deg4.yaml",
    "benchmarks/linear/systems/contraction/SOS/sos_deg8.yaml",
    "benchmarks/linear/systems/contraction/SOS/sos_deg12.yaml",
    "benchmarks/linear/systems/contraction/SOS/sos_deg24.yaml",
    "benchmarks/linear/systems/contraction/SOS/sos_deg30.yaml",
    "benchmarks/linear/systems/twotank/SOS/sos_deg4.yaml",
    "benchmarks/linear/systems/twotank/SOS/sos_deg6.yaml",
    "benchmarks/linear/systems/twotank/SOS/sos_deg8.yaml",
    "benchmarks/linear/systems/quadrotor/SOS/sos_deg4.yaml",
    "benchmarks/linear/systems/quadrotor/SOS/sos_deg6.yaml",
    "benchmarks/linear/systems/quadrotor/SOS/sos_deg8.yaml"
]

# Include your synthesis script once
include("benchmarks/barrier_synthesis.jl")

# Run all experiments in one Julia process
for yaml_file in yaml_files
    println("Running experiment: ", yaml_file)
    barrier_synthesis(yaml_file)
end

println("All SOS experiments completed!")
EOF
