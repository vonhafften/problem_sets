# Computational Economics
# Professor JF Houde
# Problem set 2
# Alex von Hafften 
# November 22, 2021

using Plots, Tables, DataFrames, CSV, StatFiles, Optim, Distributed, SharedArrays

# sets up distributed computing processes
workers()
addprocs(3)

# set working directory and reads in data 
cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_899b/ps2/")
include("data.jl");
include("model.jl");

# Initial parameter vector guess
α_0 = 0.0;
α_1 = -1.0;
α_2 = -1.0;
β = Array(fill(0.0, K_x));
γ = Array(fill(0.3, K_z));
ρ = 0.5;

include("quadrature_integration.jl");

@time likelihoods = likelihood(γ, β, ρ, α_0, α_1, α_2, t, x, z)