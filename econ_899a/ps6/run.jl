####################################################################################################################
# ECON 899A Computational Economics
# Problem Set 6

# Hopenhayn and Rogerson (1993)
# Alex von Hafften
# October 20, 2021

# This file runs the model.
####################################################################################################################

using Plots, DataFrames

include("model.jl");

####################################################################################################################
########################################## Standard ################################################################
####################################################################################################################

price = Solve_price()
results = Solve_firm_problem(Initialize(price))
M = Solve_M(results)
μ = compute_μ(results, M)
aggregate_labor = compute_labor_demand(results, M, μ)

println("Print level: ", price)
println("Mass of incumbants: ", sum((1 .- results.x) .* μ))
println("Mass of entrants: ", M)
println("Mass of entrants: ", sum(results.x .* μ))
println("Aggregate labor: ", aggregate_labor)
println("Labor of Incumbants: ", sum(results.N_d .* μ))
println("Labor of Entrants: ", M * sum(results.N_d .* ν))
println("Fraction of Labor in Entrants: ", M * sum(results.N_d .* ν)/ aggregate_labor)