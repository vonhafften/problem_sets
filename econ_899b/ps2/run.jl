# Computational Economics
# Professor JF Houde
# Problem set 2
# Alex von Hafften 
# November 22, 2021

using Plots, Tables, DataFrames, CSV, StatFiles, Optim

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_899b/ps2/")

include("data.jl");
include("integration.jl");
include("model.jl");

α_0 = 0.0
α_1 = -1.0
α_2 = -1.0
β = Array(fill(0.0, K_x))
γ = Array(fill(0.3, K_z))
ρ = 0.5

@elapsed log_likelihood(γ, β, ρ, α_0, α_1, α_2, t, x, z)