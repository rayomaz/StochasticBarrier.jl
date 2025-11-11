# using Revise, BenchmarkTools
using StochasticBarrierFunctions, LazySets
# using Mosek, MosekTools
using YAXArrays, NetCDF, MAT
using StochasticBarrierFunctions.Data

# System
number_hypercubes = 1250
σ = [0.01, 0.01, 0.01, 0.01]

filename = "../../data/unicycle/dynamics_$number_hypercubes.nc"
dataset = open_dataset(joinpath(@__DIR__, filename))
Xs = load_dynamics(dataset)

system = AdditiveGaussianUncertainPWASystem(Xs, σ)

initial_region = Hyperrectangle([-0.01, -0.01, -0.01, -0.01], [0.01, 0.01, 0.01, 0.01])
obstacle_region = EmptySet(4)

# Optimize: baseline 1 (sos)
@time res_sos = synthesize_barrier(SumOfSquaresAlgorithm(barrier_degree = 4), system, initial_region, obstacle_region)

println("Unicycle model verified.")
