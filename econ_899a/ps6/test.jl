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

results = Initialize()
P = Primitives();

@unpack s_grid, θ, β = P;

# test compute_labor_demand
results.N_d = compute_labor_demand.(1.0, s_grid, θ) # Should roughly be [1.3e-9, 10, 60, 300, 1000]

# test compute_static_profit
results.π = compute_static_profit.(1.0, s_grid, results.N_d, θ, results.c_f) # Should be increasing

# test Exit_Bellman
(results.x, results.W) = Exit_Bellman(P, results)
(results.x, results.W) = Exit_Bellman(P, results)

# test Solve_firm_problem
results = Initialize();

Solve_firm_problem(results);

results.x # a price of 1.0 is so high that no firms exit
results.π
results.W # all franchise values are positive

####################################################################################################################
############################################### Test solving for price #############################################
####################################################################################################################

results = Initialize();

compute_entry_condition(results.W, 1.0)

Solve_price(results)

####################################################################################################################
############################################### Recreate first appendix figure #####################################
####################################################################################################################

@unpack c_e = Primitives()

results = Initialize();
Solve_firm_problem(results)

price_grid = .25:.001:1.25
ec_grid = collect(copy(price_grid))

for (i, price) = enumerate(price_grid)
    results.p = price
    Solve_firm_problem(results)
    ec_grid[i] = compute_entry_condition(results.W, price) .+ c_e
end

plot(price_grid, ec_grid)
plot!(price_grid, fill(c_e, length(price_grid)))


####################################################################################################################
############################################### Test compute μ based on computing μ using t_star operator ##########
####################################################################################################################

results = Initialize()
Solve_price(results);

results.M = 2.3

compute_μ(results)
compute_μ_T_star(results)

####################################################################################################################
############################################### Test solving μ based on M ##########################################
####################################################################################################################

# Verify nonnegative over wide range of M
results = Initialize();
Solve_price(results);

M_grid = 1.0:0.1:5.0

for M in M_grid
    results.M = M
    compute_μ(results)
    if sum(results.μ .< 0) != 0  # Throws a flag if negative μ value
        error("Invalid μ")
    end
end

####################################################################################################################
############################################### Recreate second appendix figure #####################################
####################################################################################################################

M_grid = 1.0:0.1:5.0
N_s_grid = zeros(length(M_grid))
N_d_grid = zeros(length(M_grid))

results = Initialize();
Solve_price(results);

for (i, M) = enumerate(M_grid)
    results.M = M
    compute_μ(results)
    N_s_grid[i] = compute_labor_supply(results)
    N_d_grid[i] = compute_labor_demand(results)
end

plot(M_grid, N_s_grid)
plot!(M_grid, N_d_grid)

####################################################################################################################
############################################### Test solve M #######################################################
####################################################################################################################

results = Initialize();

Solve_price(results);

Solve_M(results);

results.c_f
results.p
results.N_d 
results.π
results.x
results.W
results.M
results.μ
results.Π
results.L_d
results.L_s

####################################################################################################################
############################################### Test Exit_Bellman_random ###########################################
####################################################################################################################

P = Primitives()
results = Initialize(;c_f = 10.0, α = 1.0)

Exit_Bellman_random(P, R)

Solve_firm_problem(results)

results.x