#!/bin/bash
# Usage: stochasticbarrier sos
# Single call runs all SOS and ProTECT experiments

if [ "$1" != "sos" ]; then
    echo "Usage: $0 sos"
    exit 1
fi

# -------------------------------
# Julia SOS experiments
# -------------------------------
JULIA_PROJECT="/StochasticBarrierFunctions"
echo "Starting Julia SOS experiments ..."

julia --project="$JULIA_PROJECT" <<'EOF'
# Define the YAML files for SOS
yaml_files = [

    # ---------------------Linear Systems---------------------#
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
    "benchmarks/linear/systems/quadrotor/SOS/sos_deg8.yaml",

    #---------------------Polynomial Systems---------------------#
    "benchmarks/polynomial/systems/thermostat/SOS/sos_deg2.yaml",
    "benchmarks/polynomial/systems/thermostat/SOS/sos_deg4.yaml",
    "benchmarks/polynomial/systems/thermostat/SOS/sos_deg6.yaml",
    "benchmarks/polynomial/systems/thermostat/SOS/sos_deg8.yaml",
    "benchmarks/polynomial/systems/thermostat/SOS/sos_deg12.yaml",

    "benchmarks/polynomial/systems/oscillator/SOS/sos_deg2.yaml",
    "benchmarks/polynomial/systems/oscillator/SOS/sos_deg4.yaml",
    "benchmarks/polynomial/systems/oscillator/SOS/sos_deg6.yaml",
    "benchmarks/polynomial/systems/oscillator/SOS/sos_deg8.yaml",
    "benchmarks/polynomial/systems/oscillator/SOS/sos_deg12.yaml"

    "benchmarks/polynomial/systems/room/SOS/sos_deg2.yaml",
    "benchmarks/polynomial/systems/room/SOS/sos_deg4.yaml",
    "benchmarks/polynomial/systems/room/SOS/sos_deg6.yaml",
    "benchmarks/polynomial/systems/room/SOS/sos_deg8.yaml",
    "benchmarks/polynomial/systems/room/SOS/sos_deg12.yaml"
]

# Include the barrier synthesis once
include("benchmarks/barrier_synthesis.jl")

# Run all Julia experiments
for yaml_file in yaml_files
    println("Running Julia experiment: ", yaml_file)
    barrier_synthesis(yaml_file)
end

println("All Julia SOS experiments completed!")
EOF
