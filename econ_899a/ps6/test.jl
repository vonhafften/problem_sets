####################################################################################################################
# ECON 899A Computational Economics
# Problem Set 6

# Hopenhayn and Rogerson (1993)
# Alex von Hafften
# October 20, 2021

# This file contains tests.
####################################################################################################################

include("model.jl")

####################################################################################################################
############################################### Test firm problem function #########################################
####################################################################################################################

R = Initialize(1.0);
P = Primitives(1.0);

@unpack s_grid, θ, β, c_f = P;

# test compute_labor_demand
R.N_d = compute_labor_demand.(1.0, s_grid, θ) # Should roughly be [1.3e-9, 10, 60, 300, 1000]

# test compute_static_profit
R.π = compute_static_profit.(1.0, s_grid, R.N_d, θ, c_f) # Should be increasing

# test Exit_Bellman
(R.x, R.W) = Exit_Bellman(P, R)
(R.x, R.W) = Exit_Bellman(P, R)

# test Solve_firm_problem
R = Initialize(1.0);
R = Solve_firm_problem(R);

R.x
R.π
R.W

####################################################################################################################
############################################### Test solving for price #############################################
####################################################################################################################

@elapsed compute_entry_condition(1.0)

@elapsed price = Solve_price()

Solve_firm_problem(Initialize(price))