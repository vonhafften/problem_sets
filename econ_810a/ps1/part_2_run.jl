# Kaplan and Violante (2010)
# Alex von Hafften
# January 31, 2022

# ECON 810A Advanced Macro Theory
# Problem set 1 - Part 2

# This file run analysis of the model.

# Load libraries
using Parameters, Plots

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_810a/ps1/")

####################################################################################
############### Solve model ########################################################
####################################################################################

include("part_2_model.jl")

results = Initialize();
Solve_HH_problem!(results; progress = true)

####################################################################################
############################# Simulate model #######################################
####################################################################################

include("part_2_simulation.jl")

simulation_results = Initialize_Simulation(5000)
simulate_shocks!(simulation_results)
compute_assets!(simulation_results, results)


@unpack max_a = Primitives();
sum(simulation_results.assets .== max_a)/length(simulation_results.assets)

####################################################################################
############################# Main Figures #########################################
####################################################################################

# average by age
plot(mean(simulation_results.assets, dims= 1)', legend = false, title ="Mean Wealth by Model Age");
xlabel!("Model Age");
ylabel!("Dollars")
savefig("mean_wealth.png")

# variance by age
plot(var(simulation_results.consumption, dims= 1)', legend = false, title = "Variance of Consumption by Model Age");
xlabel!("Model Age");
ylabel!("Dollars^2")
savefig("var_consumption.png")

####################################################################################
############################# Estimates ############################################
####################################################################################

# returns (var_ε, var_ζ, pass_ε, pass_ζ)
bbp_coefficient = compute_bbp_coefficients(simulation_results)

####################################################################################
############################# Other Figures ########################################
####################################################################################

# Check out value and policy functions
plot(results.value_function[35, :, 1, 1]);
plot!(results.value_function[35, :, 5, 5])

plot(results.value_function[5, :, 1, 1]);
plot!(results.value_function[5, :, 5, 5])

plot(results.value_function[35, :, 3, 3]);
plot!(results.value_function[34, :, 3, 3])

# plot(results.policy_function[5, :, 1, 1]);
# plot!(results.policy_function[5, :, 5, 5])

# plot(results.policy_function[10, :, 1, 1]);
# plot!(results.policy_function[10, :, 5, 5])

plot(results.policy_function[30, :, 1, 1]);
plot!(results.policy_function[30, :, 5, 5])


plot(results.policy_function[35, :, 3, 3]);
plot!(results.policy_function[34, :, 3, 3]);
plot!(results.policy_function[33, :, 3, 3])


plot(mean(simulation_results.income, dims= 1)', legend = false)
plot(mean(simulation_results.consumption, dims= 1)', legend = false)

plot(var(simulation_results.assets, dims= 1)', legend = false)
plot(var(simulation_results.income, dims= 1)', legend = false)