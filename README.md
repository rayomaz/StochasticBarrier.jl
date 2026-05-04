# StochasticBarrier.jl

We present *StochasticBarrier.jl*, an open-source Julia-based toolbox for generating Stochastic Barrier Functions (SBFs) for safety verification of discrete-time stochastic systems with additive Gaussian noise. The tool supports linear, polynomial, and piecewise affine (PWA) uncertain dynamics. The toolbox implements a Sum-of-Squares (SOS) optimization approach, as well as methods based on piecewise constant (PWC) functions. For the class of PWC-SBFs, three engines are offered based on: (1) DUAL Linear Programming, (2) Counter Example Guided (CEGS) Linear Programming, and (3) Projected Gradient Descent (GD).

## Purpose of this code
This code generates results of the benchmarks presented in Table (1) and Table (2) of the toolbox paper.
A total of six experiments are included for benchmarking SOS:
1.  **Contraction Map**
2.  **Two Tank**
3.  **Quadrotor**
4.  **Thermostat**
5.  **Oscillator**
6.  **Room Temperature**

A total of three experiments are included for benchmarking PWC:
1.  **Contraction Map**
2.  **Pendulum**
3.  **Unicycle**

## Repeat Experiments
| **`Linux`** | **`Mac OS X`** | **`Windows`** |
|-----------------|---------------------|-------------------------|

Read the description below for repeatability of all the experiments.

### Mosek License
All SOS experiments use Mosek as the SDP solver. A valid Mosek license file is therefore required.

# Obtain your license
Obtain `mosek.lic` from https://www.mosek.com/; Mosek provides free academic licenses. 

# Place the license
Put the file at
- `StochasticBarrierFunctions/benchmarks/mosek/mosek.lic`

The license is bind-mounted into the docker container at runtime (see the `docker run` command below).

### Docker Image
A pre-built image tarball is distributed alongside this repository. Load it with
```sh
sudo docker load -i stochastic_barrier.tar
```

If you prefer to build it yourself the `Dockerfile` at the repository root is still available: `sudo docker build -t stochastic_barrier .`.

To start a container that writes benchmark CSVs back to the host, bind-mount both the Mosek license and a results directory:

```sh
mkdir -p results
sudo docker run --rm -it --name StochasticBarrier \
    -v "$(pwd)/StochasticBarrierFunctions/benchmarks/mosek/mosek.lic:/mosek/mosek.lic:ro" \
    -v "$(pwd)/results:/results" \
    stochastic_barrier
```

Inside the container `MOSEKLM_LICENSE_FILE=/mosek/mosek.lic` and `SB_RESULTS_DIR=/results` are pre-set; the benchmarks will append CSV rows to `/results/julia_results.csv`.

## Run through bash

Inside the container, use:

```sh
stochasticbarrier sos                  # To run the SOS barrier benchmark
stochasticbarrier pwc                  # To run the PWC barrier benchmark
```

## Parse results into Table 3 / Table 4

After running StochasticBarrier.jl (Julia), PRoTECT (Python) and StochasticBarrierFunc. (MATLAB) benchmarks, run the parser on the host to produce the paper's two summary tables:

```sh
julia --project=StochasticBarrierFunctions/benchmarks \
    StochasticBarrierFunctions/benchmarks/parse_results.jl results
```
The above command requires Julia with only stdlibs installed.

The command reads `results/{julia,protect,sostools}_results.csv` and writes `results/table3.csv` and `results/table4.csv` (wide format matching the tables in the paper; missing rows are reported as `OM`). 