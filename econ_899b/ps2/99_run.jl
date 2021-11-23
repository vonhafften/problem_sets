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

# loads scripts
@everywhere include("01_toolbox.jl");
include("02_likelihood.jl");

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
u_0_p, u_1_p, u_2_p = initialize_ghk(;use_halton = false)
ε_0_p, ε_1_p, ε_2_p = initialize_accept_reject(ρ; use_halton = false)

# Part 1 - Quadrature Method Integration
@time likelihoods_quadrature = likelihood(α_0, α_1, α_2, β, γ, ρ, t, x, z, q_grids[1], q_grids[2], u_0_h, u_1_h, u_2_h, ε_0_h, ε_1_h, ε_2_h; method = "quadrature")

# Part 2 - GHK Method
@time likelihoods_ghk        = likelihood(α_0, α_1, α_2, β, γ, ρ, t, x, z, q_grids[1], q_grids[2], u_0_h, u_1_h, u_2_h, ε_0_h, ε_1_h, ε_2_h; method = "ghk")
@time likelihoods_ghk_pseudo = likelihood(α_0, α_1, α_2, β, γ, ρ, t, x, z, q_grids[1], q_grids[2], u_0_p, u_1_p, u_2_p, ε_0_p, ε_1_p, ε_2_p; method = "ghk")

# Part 3 - Accept-Reject Method
@time likelihoods_accept_reject        = likelihood(α_0, α_1, α_2, β, γ, ρ, t, x, z, q_grids[1], q_grids[2], u_0_h, u_1_h, u_2_h, ε_0_h, ε_1_h, ε_2_h; method = "accept_reject")
@time likelihoods_accept_reject_pseudo = likelihood(α_0, α_1, α_2, β, γ, ρ, t, x, z, q_grids[1], q_grids[2], u_0_p, u_1_p, u_2_p, ε_0_p, ε_1_p, ε_2_p; method = "accept_reject")

# save results
p4_result = copy(df)

p4_result[!, :likelihood_quadrature]           = likelihoods_quadrature
p4_result[!, :likelihood_ghk]                  = likelihoods_ghk
p4_result[!, :likelihood_ghk_pseudo]           = likelihoods_ghk_pseudo
p4_result[!, :likelihood_accept_reject]        = likelihoods_accept_reject
p4_result[!, :likelihood_accept_reject_pseudo] = likelihoods_accept_reject_pseudo

CSV.write("p4_result.csv", p4_result)

# Part 5 - LBFGS optimization
# My optimization code did not converge here are the parameters with the highest log-likelihood so far...
# Initial parameter vector guess
# Log-likelihood: -11480.20991042999
α_0 = -2.678674561062139
α_1 = -1.2476328255418148
α_2 = -1.0356907807723135
β = [-0.149411873515521, 0.2136948477122363, 0.3937060642915554, 0.2425304814585507, 0.06736233478802248, -0.014366052643485401, -0.0010864692187960538, 0.030252205164349295, 0.18166591312677577, 0.17140386077054118, 0.21291304329070787, 0.4528297037610101, 0.14078126937653634, -0.0519203088138539, -0.03813958833558798]
γ = -0.021475182529755162
ρ = 0.12926220640818922

optimization_results = optimize(θ -> -log_likelihood(θ, t, x, z, q_grids[1], q_grids[2], u_0_h, u_1_h, u_2_h, ε_0_h, ε_1_h, ε_2_h), 
                                vcat(α_0, α_1, α_2, β, γ, ρ),
                                LBFGS())
