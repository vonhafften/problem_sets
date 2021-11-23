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
# Initial parameter vector guess
# Log-likelihood: -11536.256813515432
α_0 = -0.6681692422928824
α_1 = -0.4430814561230035
α_2 = -0.2742889000441864
β = [-0.16589035210469164, 0.16438566676493688, 0.2727737054741467, 0.1944853339742954, 0.023948722238739402, 0.15712086933455632, -0.35010978938640397, -0.014863922928109171, 0.11657324499118413, 0.10197009559120476, 0.17976736557309758, 0.3703296969761364, 0.13241512183563, 0.021063897032843865, -0.009041343050043449]
γ = -0.03800314475586985
ρ = -0.159223672078737

optimization_results = optimize(θ -> -log_likelihood(θ, t, x, z, q_grids[1], q_grids[2], u_0_h, u_1_h, u_2_h, ε_0_h, ε_1_h, ε_2_h), 
                                vcat(α_0, α_1, α_2, β, γ, ρ),
                                LBFGS())
