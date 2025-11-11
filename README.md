# StochasticBarrier.jl

We present *StochasticBarrier.jl*, an open-source Julia-based toolbox for generating Stochastic Barrier Functions (SBFs) for safety verification of discrete-time stochastic systems with additive Gaussian noise. The tool supports linear, polynomial, and piecewise affine (PWA) uncertain dynamics. The toolbox implements a Sum-of-Squares (SOS) optimization approach, as well as  methods based on piecewise constant (PWC) functions. For the class of of PWC-SBFs, three engines are offered based on: (1) DUAL Linear Programming, (2) Counter Example Guided (CEGS) Linear Programming, and (3) Projected Gradient Descent (GD).

## Purpose of this code
This code generates results of the benchmarks presented in Table (1) and Table (2) of the toolbox paper.
A total of six experiments are included for benchmarking SOS:
1.  **Contraction Map**
2.  **Two Tank**
3.  **Quadrotor**
4.  **Thermostat**
5.  **Oscillator**
6.  **Room Temperature**

A total of four experiments are included for benchmarking PWC:
1.  **Contraction Map**
2.  **Pendulum**
3.  **Unicycle**

## Repeat Experiments
| **`Linux`** | **`Mac OS X`** | **`Windows`** |
|-----------------|---------------------|-------------------------|

Read the description below for repeatability of all the experiments.

### Docker Image
The Dockerfile is provided in the main folder. Build this docker file to obtain all the required Julia packages, as specified in the Project.toml. To build the docker image, navigate into the main folder and run the following command 
```sh
sudo docker build -t stochastic_barrier .
```

To start a container 

```sh
sudo docker run -it --name StochasticBarrier stochastic_barrier
```

## Run through bash

Use the following commands to run the benchmarks.

```sh
stochasticbarrier sos                  # To run the SOS barrier benchmark
stochasticbarrier pwc                  # To run the PWC barrier benchmark
```
