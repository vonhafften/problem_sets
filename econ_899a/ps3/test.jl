# Conesa and Krueger (1999)
# Alex von Hafften
# September 22, 2021

# ECON 899A Computational Economics
# Problem set 3

# This file calls model.jl and is used for testing each function.

using Plots

################################################################################
########################## Tests Solve_retiree_problem() ####################
################################################################################

include("model.jl")

@unpack a_grid = Primitives()

results = Initialize()

Solve_retiree_problem(results; progress = true)

plot(a_grid, results.value_function[66, :, 1])
plot(a_grid, results.value_function[47, :, 1])
plot(a_grid, results.value_function[60:66, :, 1]')
plot(a_grid, results.value_function[50:59, :, 1]')
plot(a_grid, results.value_function[47:49, :, 1]')

plot(a_grid, results.policy_function[66, :, 1])
plot(a_grid, results.policy_function[65, :, 1])
plot(a_grid, results.policy_function[47, :, 1])
plot(a_grid, results.policy_function[60:66, :, 1]')
plot(a_grid, results.policy_function[50:59, :, 1]')
plot(a_grid, results.policy_function[47:49, :, 1]')
