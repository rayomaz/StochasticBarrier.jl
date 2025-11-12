using Revise, BenchmarkTools
using StochasticBarrierFunctions, LazySets, LinearAlgebra
using YAXArrays, NetCDF

# System
F = 0.50
dim = 2

A = F*I(dim)
b = zeros(dim)
σ = 0.1*ones(dim)

state_space = Hyperrectangle(low=-1.0*ones(dim), high=0.5*ones(dim))

system = AdditiveGaussianLinearSystem(A, b, σ, state_space)

# Initial range and obstacle space
initial_region = Hyperrectangle(low=[-0.05, -0.05], high=[0.05, 0.05])

obstacle_region = EmptySet(dim)

# Set horizon
N = 10

# Optimize: baseline 1 (sos)
@time res_sos = synthesize_barrier(SumOfSquaresAlgorithm(barrier_degree = 30), system, initial_region, obstacle_region; time_horizon = N)

println("Single mode model verified.")
