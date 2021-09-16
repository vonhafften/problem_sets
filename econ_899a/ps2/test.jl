# Huggett (1993)
# Alex von Hafften
# September 16, 2021

# ECON 899A Computational Economics
# Problem set 2

# This file calls model.jl and is used for testing each component.

include("model.jl")

################################################################################
########################## Tests HH_Bellman() ##################################
################################################################################

@unpack a_grid = Primitives()

results = Initialize()

v_next = HH_Bellman(results)

plot(a_grid, results.value_function) # Should be lines at zero.
plot(a_grid, results.policy_function) # Should look like noise due to being zoomed in so much

results.value_function = v_next

v_next = HH_Bellman(results)

plot(a_grid, results.value_function) # Should look more like the value function we expect
plot(a_grid, results.policy_function) # Same as above but for policy function

################################################################################
########################## Tests Solve_HH_problem() ############################
################################################################################

@unpack a_grid = Primitives()

results = Initialize()

Solve_HH_problem(results)

plot(a_grid, results.value_function)
plot(a_grid, results.policy_function)

################################################################################
########################## Tests Solve_HH_problem() ############################
################################################################################
