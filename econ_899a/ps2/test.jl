# Huggett (1993)
# Alex von Hafften
# September 16, 2021

# ECON 899A Computational Economics
# Problem set 2

# This file calls model.jl and is used for testing each function.

using Plots

################################################################################
########################## Tests HH_Bellman() ##################################
################################################################################

include("model.jl")

@unpack a_grid = Primitives()

results = Initialize()

v_next = HH_Bellman(results; progress = true)

plot(a_grid, results.value_function) # Should be lines at zero.
plot(a_grid, results.policy_function) # Should look like noise due to being zoomed in so much

results.value_function = v_next

v_next = HH_Bellman(results; progress = true)

plot(a_grid, results.value_function) # Should look more like the value function we expect
plot(a_grid, results.policy_function) # Same as above but for policy function

################################################################################
########################## Tests Solve_HH_problem() ############################
################################################################################

include("model.jl")

@unpack a_grid = Primitives()

results = Initialize()

Solve_HH_problem(results; progress = true)

plot(a_grid, results.value_function)
plot(a_grid, results.policy_function)

################################################################################
########################## Tests Apply_policy_function_to_μ() ##################
################################################################################

include("model.jl")

@unpack a_grid = Primitives()

results = Initialize()

Solve_HH_problem(results)

μ_next = Apply_policy_function_to_μ(results; progress = true)

plot(a_grid, results.μ)
plot(a_grid, μ_next)

################################################################################
########################## Tests Solve_invariant_μ() ###########################
################################################################################

include("model.jl")

@unpack a_grid = Primitives()

results = Initialize()

Solve_HH_problem(results)

Solve_invariant_μ(results; progress = true)

plot(a_grid, results.μ)

################################################################################
########################## Tests update_price() ################################
################################################################################

include("model.jl")

@unpack a_grid = Primitives()

results = Initialize()

Solve_HH_problem(results)

Solve_invariant_μ(results)

println(results.q)

Update_price(results)

println(results.q)

################################################################################
########################## Tests Solve_model() #################################
################################################################################

include("model.jl")

@unpack a_grid = Primitives()

results = Initialize()

Solve_model(results)

plot(a_grid, results.value_function, labels = ["Employed" "Unemployed"])
plot(a_grid, [results.policy_function a_grid], labels = ["Employed" "Unemployed" "45° Line"])
plot(a_grid, results.μ, labels = ["Employed" "Unemployed"])
plot(a_grid, cumsum(results.μ; dims = 1), labels = ["Employed" "Unemployed"])
