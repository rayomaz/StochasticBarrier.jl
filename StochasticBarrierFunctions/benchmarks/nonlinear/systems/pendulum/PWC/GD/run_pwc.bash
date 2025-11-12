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
# Define the YAML files for SOS and PWC
yaml_files = [

    # ---------------------SOS Linear System---------------------#
    # "benchmarks/linear/systems/contraction2/SOS/sos_deg4.yaml",
    # "benchmarks/linear/systems/contraction2/SOS/sos_deg8.yaml",
    # "benchmarks/linear/systems/contraction2/SOS/sos_deg20.yaml",
    # "benchmarks/linear/systems/contraction2/SOS/sos_deg30.yaml",

    # ---------------------PWC Linear System---------------------#
    #"benchmarks/linear/systems/contraction2/PWC/DUAL/contraction_dual_0.20.yaml",
    #"benchmarks/linear/systems/contraction2/PWC/DUAL/contraction_dual_0.10.yaml",
    #"benchmarks/linear/systems/contraction2/PWC/DUAL/contraction_dual_0.09.yaml",
    #"benchmarks/linear/systems/contraction2/PWC/DUAL/contraction_dual_0.08.yaml",

    #"benchmarks/linear/systems/contraction2/PWC/CEGIS/contraction_cegis_0.20.yaml",
    #"benchmarks/linear/systems/contraction2/PWC/CEGIS/contraction_cegis_0.10.yaml",
    #"benchmarks/linear/systems/contraction2/PWC/CEGIS/contraction_cegis_0.09.yaml",
    #"benchmarks/linear/systems/contraction2/PWC/CEGIS/contraction_cegis_0.08.yaml",

    #"benchmarks/linear/systems/contraction2/PWC/GD/contraction_gd_0.20.yaml",
    #"benchmarks/linear/systems/contraction2/PWC/GD/contraction_gd_0.10.yaml",
    #"benchmarks/linear/systems/contraction2/PWC/GD/contraction_gd_0.09.yaml",
    #"benchmarks/linear/systems/contraction2/PWC/GD/contraction_gd_0.08.yaml",

    # ---------------------SOS NNDM System---------------------#
    # "benchmarks/nonlinear/systems/pendulum/SOS/sos_120.yaml",
    # "benchmarks/nonlinear/systems/pendulum/SOS/sos_240.yaml",
    # "benchmarks/nonlinear/systems/pendulum/SOS/sos_480.yaml,

    # ---------------------PWC NNDM System---------------------#
    "benchmarks/nonlinear/systems/pendulum/PWC/DUAL/pwc_dual_120.yaml",
    "benchmarks/nonlinear/systems/pendulum/PWC/DUAL/pwc_dual_240.yaml",
    "benchmarks/nonlinear/systems/pendulum/PWC/DUAL/pwc_dual_480.yaml"

    # ---------------------SOS Hybrid System---------------------#

    # ---------------------PWC Hybrid System---------------------#


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
