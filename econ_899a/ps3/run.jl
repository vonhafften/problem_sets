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

include("model.jl")

@unpack a_grid = Primitives()

results = Initialize()

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
################################## Exercise 2 ##################################
################################################################################
