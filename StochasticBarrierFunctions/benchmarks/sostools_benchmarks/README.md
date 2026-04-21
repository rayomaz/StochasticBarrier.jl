Read the description below for repeatability of the SOSTOOLS experiments.

### MATLAB Setup
These experiments require **MATLAB** with **Symbolic Math Toolbox** to be installed on your local machine.

Dependencies **SOSTOOLS** and **MOSEK** will be installed through the provided bash script.

### Output directory
Each experiment appends a row to `$SB_RESULTS_DIR/sostools_results.csv`. If `SB_RESULTS_DIR` is not set the script falls back to
`StochasticBarrierFunctions/benchmarks/sostools_benchmarks/results/`.

Point it at the same directory used for the StochasticBarrier.jl and PRoTECT runs so the results parser can pick them all up:

```sh
export SB_RESULTS_DIR=$(pwd)/../../../results   # or an absolute path
```

## Run through bash

Make the bash file executable:
```sh
chmod +x sostoolsbench.sh
```

Then run:

```sh
./sostoolsbench.sh
```

Each experiment appends a row to `/results/sostools_results.csv`. See README.md in the repo root for how to combine with the results from StochasticBarrier.jl and StochasticBarrierFunc. for a combined table.