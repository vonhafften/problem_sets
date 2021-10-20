####################################################################################################################
# ECON 899A Computational Economics
# Problem Set 6

# Hopenhayn and Rogerson (1993)
# Alex von Hafften
# October 20, 2021

# This file contains tests.
####################################################################################################################

using Plots

include("model.jl");

####################################################################################################################
############################################### Test firm problem function #########################################
####################################################################################################################

results = Initialize(1.0);
P = Primitives();

@unpack s_grid, θ, β, c_f = P;

# test compute_labor_demand
results.N_d = compute_labor_demand.(1.0, s_grid, θ) # Should roughly be [1.3e-9, 10, 60, 300, 1000]

# test compute_static_profit
results.π = compute_static_profit.(1.0, s_grid, results.N_d, θ, c_f) # Should be increasing

# test Exit_Bellman
(results.x, results.W) = Exit_Bellman(P, results)
(results.x, results.W) = Exit_Bellman(P, results)

# test Solve_firm_problem
results = Solve_firm_problem(Initialize(1.0));

results.x # a price of 1.0 is so high that no firms exit
results.π
results.W # all franchise values are positive

####################################################################################################################
############################################### Test solving for price #############################################
####################################################################################################################

compute_entry_condition(1.0)

@elapsed price = Solve_price()

results = Solve_firm_problem(Initialize(price))

####################################################################################################################
############################################### Recreate first appendix figure #####################################
####################################################################################################################

@unpack c_e = Primitives()

price_grid = .25:.001:1.25
ec_grid = compute_entry_condition.(price_grid) .+ c_e

plot(price_grid, ec_grid)
plot!(price_grid, fill(c_e, length(price_grid)))


####################################################################################################################
############################################### Test compute μ based on computing μ using t_star operator ##########
####################################################################################################################

results = Solve_firm_problem(Initialize(Solve_price()))

compute_μ(results, 2.3)
compute_μ_T_star(results, 2.3)

####################################################################################################################
############################################### Test solving μ based on M ##########################################
####################################################################################################################

# Verify nonnegative over wide range of M
M_grid = 1.0:0.1:5.0

results = Solve_firm_problem(Initialize(Solve_price()))

for M in M_grid
    μ = compute_μ(results, M)
    if sum(μ .< 0) != 0  # Throws a flag if negative μ value
        error("Invalid μ")
    end
end

####################################################################################################################
############################################### Recreate second appendix figure #####################################
####################################################################################################################

M_grid = 1.0:0.1:5.0
N_s_grid = zeros(length(M_grid))
N_d_grid = zeros(length(M_grid))

results = Solve_firm_problem(Initialize(Solve_price()))

for (i, M) = enumerate(M_grid)
    μ = compute_μ(results, M)
    N_s_grid[i] = compute_labor_supply(results, M, μ)
    N_d_grid[i] = compute_labor_demand(results, M, μ)
end

plot(M_grid, N_s_grid)
plot!(M_grid, N_d_grid)

####################################################################################################################
############################################### Text solve M #######################################################
####################################################################################################################

results = Solve_firm_problem(Initialize(Solve_price()))

Solve_M(results)

@elapsed Solve_M(Solve_firm_problem(Initialize(Solve_price())))