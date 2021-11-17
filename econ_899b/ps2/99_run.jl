# Computational Economics
# Professor JF Houde
# Problem set 2
# Alex von Hafften 
# November 22, 2021

using Plots, Optim, Distributed, SharedArrays

# sets up distributed computing processes
procs() # should be 1
addprocs(3)
procs() # should be 1:4

# set working directory and reads in data 
@everywhere cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_899b/ps2/")
include("data.jl");
include("model.jl");

# Initial parameter vector guess
α_0 = 0.0;
α_1 = -1.0;
α_2 = -1.0;
β = Array(fill(0.0, K_x));
γ = Array(fill(0.3, K_z));
ρ = 0.5;

# Part 1 - Quadrature Method Integration
@everywhere include("quadrature_integration.jl");
@elapsed quadrature_ll = log_likelihood(γ, β, ρ, α_0, α_1, α_2, t, x, z)

# Part 2 - GHK Method Integration

# Part 3 - Accept-Reject Method Integration

# Part 5 - Gradiant-free optimization
# optimize(θ -> -log_likelihood(θ, t, x, z), vcat(γ, β, ρ, α_0, α_1, α_2))