# Computational Economics
# Professor JF Houde
# Problem set 2
# Alex von Hafften 
# November 22, 2021

using Plots, Optim, Distributed, SharedArrays

# sets up distributed computing processes
workers()
addprocs(3)

# loads other scripts
include("01_data.jl");
@everywhere include("02_toolbox.jl")
include("03_likelihood.jl");

# Initial parameter vector guess
α_0 =  0.0;
α_1 = -1.0;
α_2 = -1.0;
β = Array(fill(0.0, K_x));
γ = Array(fill(0.3, K_z));
ρ = 0.5;

@time likelihoods = likelihood(γ, β, ρ, α_0, α_1, α_2, t, x, z; method = "quadrature")

@time likelihoods = likelihood(γ, β, ρ, α_0, α_1, α_2, t, x, z; method = "accept_reject")