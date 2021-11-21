# Computational Economics
# Professor JF Houde
# Problem set 2
# Alex von Hafften 
# November 22, 2021

using Plots, Optim, DataFrames, CSV, Tables, Distributed, SharedArrays, ProgressMeter, Statistics, StatFiles, Random

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_899b/ps2/")

# sets up distributed computing processes
workers()
addprocs(3)

# load data
df = DataFrame(load("PS2/Mortgage_performance_data.dta"))

# Define x, y, and z.
df_x = select(df, 
              :score_0, :rate_spread, :i_large_loan, :i_medium_loan, 
              :i_refinance, :age_r, :cltv, :dti, :cu,  :first_mort_r, :i_FHA,
              :i_open_year2, :i_open_year3, :i_open_year4, :i_open_year5)
df_t = select(df, :duration)
df_z = select(df, :score_0, :score_1, :score_2)    

# Define matrices.
x = Float64.(Array(df_x))
z = Float64.(Array(df_z))
t = Float64.(Array(df_t))

# Initial parameter vector guess
α_0 =  0.0;
α_1 = -1.0;
α_2 = -1.0;
β   = Array(fill(0.0, size(x)[2]));
γ   = Array(fill(0.3, size(z)[2]));
ρ   = 0.5;

# loads scripts
@everywhere include("01_toolbox.jl");
include("02_likelihood.jl");

# Part 1 - Quadrature Method Integration
@time likelihoods_quadrature = likelihood(γ, β, ρ, α_0, α_1, α_2, t, x, z; method = "quadrature")

# Part 2 - GHK Method
@time likelihoods_ghk = likelihood(γ, β, ρ, α_0, α_1, α_2, t, x, z; method = "ghk")
@time likelihoods_ghk_pseudo = likelihood(γ, β, ρ, α_0, α_1, α_2, t, x, z; method = "ghk", use_halton = false)

# Part 3 - Accept-Reject Method
@time likelihoods_accept_reject = likelihood(γ, β, ρ, α_0, α_1, α_2, t, x, z; method = "accept_reject")
@time likelihoods_accept_reject_pseudo = likelihood(γ, β, ρ, α_0, α_1, α_2, t, x, z; method = "accept_reject", use_halton = false)

# save results
p4_result = copy(df)

p4_result[!, :likelihood_quadrature] = likelihoods_quadrature
p4_result[!, :likelihoods_ghk] = likelihoods_ghk
p4_result[!, :likelihoods_ghk_pseudo] = likelihoods_ghk_pseudo
p4_result[!, :likelihoods_accept_reject] = likelihoods_accept_reject
p4_result[!, :likelihoods_accept_reject_pseudo] = likelihoods_accept_reject_pseudo

CSV.write("p4_result.csv", p4_result)

# Part 5 - Gradiant-free optimization
# optimization_results = optimize(θ -> -log_likelihood(θ, t, x, z), vcat(γ, β, ρ, α_0, α_1, α_2))