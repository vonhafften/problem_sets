# Computational Economics
# Professor JF Houde
# Problem set 2
# Alex von Hafften 
# November 22, 2021

using Plots, Optim, DataFrames, CSV, Tables, Distributed, SharedArrays, ProgressMeter, Statistics, StatFiles, Random, Distributions

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
γ   = 0.3;
ρ   = 0.5;

# quadrature points
q_grids = initialize_quadrature_integration()

# random shocks
u_0_h, u_1_h, u_2_h = initialize_ghk(;use_halton = true)
ε_0_h, ε_1_h, ε_2_h = initialize_accept_reject(ρ; use_halton = true)

# loads scripts
@everywhere include("01_toolbox.jl");
include("02_likelihood.jl");

# Part 1 - Quadrature Method Integration

@time likelihoods_quadrature_1 = likelihood(α_0, α_1, α_2, β, γ, ρ, t./t*1.0, x, z, q_grids[1], q_grids[2], u_0_h, u_1_h, u_2_h, ε_0_h, ε_1_h, ε_2_h; method = "quadrature")
@time likelihoods_quadrature_2 = likelihood(α_0, α_1, α_2, β, γ, ρ, t./t*2.0, x, z, q_grids[1], q_grids[2], u_0_h, u_1_h, u_2_h, ε_0_h, ε_1_h, ε_2_h; method = "quadrature")
@time likelihoods_quadrature_3 = likelihood(α_0, α_1, α_2, β, γ, ρ, t./t*3.0, x, z, q_grids[1], q_grids[2], u_0_h, u_1_h, u_2_h, ε_0_h, ε_1_h, ε_2_h; method = "quadrature")
@time likelihoods_quadrature_4 = likelihood(α_0, α_1, α_2, β, γ, ρ, t./t*4.0, x, z, q_grids[1], q_grids[2], u_0_h, u_1_h, u_2_h, ε_0_h, ε_1_h, ε_2_h; method = "quadrature")

likelihoods_quadrature_1 .+ likelihoods_quadrature_2 .+ likelihoods_quadrature_3 .+ likelihoods_quadrature_4

# Part 2 - GHK Method

@time likelihoods_ghk_1 = likelihood(α_0, α_1, α_2, β, γ, ρ, t./t*1.0, x, z, q_grids[1], q_grids[2], u_0_h, u_1_h, u_2_h, ε_0_h, ε_1_h, ε_2_h; method = "ghk")
@time likelihoods_ghk_2 = likelihood(α_0, α_1, α_2, β, γ, ρ, t./t*2.0, x, z, q_grids[1], q_grids[2], u_0_h, u_1_h, u_2_h, ε_0_h, ε_1_h, ε_2_h; method = "ghk")
@time likelihoods_ghk_3 = likelihood(α_0, α_1, α_2, β, γ, ρ, t./t*3.0, x, z, q_grids[1], q_grids[2], u_0_h, u_1_h, u_2_h, ε_0_h, ε_1_h, ε_2_h; method = "ghk")
@time likelihoods_ghk_4 = likelihood(α_0, α_1, α_2, β, γ, ρ, t./t*4.0, x, z, q_grids[1], q_grids[2], u_0_h, u_1_h, u_2_h, ε_0_h, ε_1_h, ε_2_h; method = "ghk")

likelihoods_ghk_1 .+ likelihoods_ghk_2 .+ likelihoods_ghk_3 .+ likelihoods_ghk_4

# Part 3 - Accept-Reject Method

@time likelihoods_accept_reject_1 = likelihood(α_0, α_1, α_2, β, γ, ρ, t./t*1.0, x, z, q_grids[1], q_grids[2], u_0_h, u_1_h, u_2_h, ε_0_h, ε_1_h, ε_2_h; method = "accept_reject")
@time likelihoods_accept_reject_2 = likelihood(α_0, α_1, α_2, β, γ, ρ, t./t*2.0, x, z, q_grids[1], q_grids[2], u_0_h, u_1_h, u_2_h, ε_0_h, ε_1_h, ε_2_h; method = "accept_reject")
@time likelihoods_accept_reject_3 = likelihood(α_0, α_1, α_2, β, γ, ρ, t./t*3.0, x, z, q_grids[1], q_grids[2], u_0_h, u_1_h, u_2_h, ε_0_h, ε_1_h, ε_2_h; method = "accept_reject")
@time likelihoods_accept_reject_4 = likelihood(α_0, α_1, α_2, β, γ, ρ, t./t*4.0, x, z, q_grids[1], q_grids[2], u_0_h, u_1_h, u_2_h, ε_0_h, ε_1_h, ε_2_h; method = "accept_reject")

likelihoods_accept_reject_1 .+ likelihoods_accept_reject_2 .+ likelihoods_accept_reject_3 .+ likelihoods_accept_reject_4
