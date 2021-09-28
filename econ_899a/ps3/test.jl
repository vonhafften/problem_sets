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

@elapsed Solve_retiree_problem(results; progress = true)

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

################################################################################
########################## Tests Solve_retiree_problem() #######################
################################################################################

include("model.jl")

@unpack a_grid = Primitives()

results = Initialize()

@elapsed Solve_retiree_problem(results)

@elapsed Solve_worker_problem(results; progress = true)

plot(a_grid, results.value_function[28, :, 2])
plot(a_grid, results.value_function[20:25, :, 1]')
plot(a_grid, results.value_function[25:30, :, 2]')
plot([a_grid a_grid], results.value_function[20, :, 1:2])
plot([a_grid a_grid], results.policy_function[41, :, 1:2])
plot(a_grid, results.labor_supply[40:45, :, 1]')

################################################################################
########################## Tests Solve_HH_problem() ############################
################################################################################

include("model.jl")

@unpack a_grid = Primitives()

results = Initialize()

@elapsed Solve_HH_problem(results)
