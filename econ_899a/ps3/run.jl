# Conesa and Krueger (1999)
# Alex von Hafften
# September 22, 2021

# ECON 899A Computational Economics
# Problem set 3

# This file calls model.jl and runs the computations

include("model.jl");

using Plots

################################################################################
################################## Exercise 1 ##################################
################################################################################

@unpack a_grid = Primitives()

results = Initialize(0.11, [3.0, 0.5], 0.42)

@elapsed Solve_HH_problem(results)

plot(a_grid,
     results.value_function[50, :, 1],
     labels = "",
     legend=:bottomright)

savefig("value_function_50.png")

plot([a_grid a_grid a_grid],
     [results.policy_function[20, :, :] a_grid],
     labels = ["High" "Low" "45° Line"],
     legend=:bottomright)

savefig("policy_function_20.png")

################################################################################
################################## Exercise 3 ##################################
################################################################################

# Benchmark
@elapsed bm_ss = Solve_model()  # converges in ~9 iterations
@elapsed bm_no_ss = Solve_model(θ = 0.0)  # converges in ~11 iterations

# No productivity risk
@elapsed riskless_ss = Solve_model(z = [0.5, 0.5]) # converges in ~12 iterations
@elapsed riskless_no_ss = Solve_model(θ = 0.0, z = [0.5, 0.5], λ = 0.1)  # converges in ~52 iterations

# Inelastic labor supply
@elapsed inelastic_l_ss = Solve_model(γ = 1.0, λ = 0.8) # converges in ~6 iterations
@elapsed inelastic_l_no_ss = Solve_model(θ = 0.0, γ = 1.0, λ = 0.8) # converges in ~7 iterations

table_1 =create_table([bm_ss, bm_no_ss,
                        riskless_ss, riskless_no_ss,
                        inelastic_l_ss, inelastic_l_no_ss])

CSV.write("table_1.csv", table_1)
