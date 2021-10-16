####################################################################################################################
# Krusell and Smith (1998)
# Alex von Hafften
# October 13, 2021

# ECON 899A Computational Economics
# Problem set 5

# This file contains test functions.
####################################################################################################################

using Plots

####################################################################################################################
############################# test simulate_Z() ####################################################################
####################################################################################################################

include("02_initialize.jl");

@unpack Π_z_star = Shocks();
@unpack T = Simulation();

n_trials = 10000
recession_share = zeros(n_trials)

for i = 1:n_trials
    Z = simulate_Z()
    recession_share[i] = count(Z .== 2)/T
end

println(Π_z_star)

histogram(recession_share) # should be centered at the stationary distribution

sum(recession_share)/n_trials

####################################################################################################################
############################# test simulate_E(Z) ###################################################################
####################################################################################################################

include("02_initialize.jl");

@unpack N, T, burn_in = Simulation();
@unpack Π_ε_star = Shocks();

n_trials = 100;

boom_employment = zeros(n_trials);
boom_unemployment = zeros(n_trials);
recession_employment = zeros(n_trials);
recession_unemployment = zeros(n_trials);

# takes a bit...
for i = 1:n_trials
    println(i)

    Z = simulate_Z()
    E = simulate_E(Z)

    E = E[:, burn_in:T]
    Z = Z[burn_in:T]

    E_boom      = E[:, Z .== 1]
    E_recession = E[:, Z .== 2]

    boom_employment[i] = count(E_boom .== 1)/(N*(T-burn_in+1))
    boom_unemployment[i] = count(E_boom .== 2)/(N*(T-burn_in+1))
    recession_employment[i] = count(E_recession .== 1)/(N*(T-burn_in+1))
    recession_unemployment[i] = count(E_recession .== 2)/(N*(T-burn_in+1))
end

println(Π_ε_star)

histogram(boom_employment)
histogram(boom_unemployment)
histogram(recession_employment)
histogram(recession_unemployment)

sum(boom_employment)/n_trials
sum(boom_unemployment)/n_trials
sum(recession_employment)/n_trials
sum(recession_unemployment)/n_trials

####################################################################################################################
############################# test Bellman(R) ######################################################################
####################################################################################################################

include("03_solve_value_function.jl");

@unpack k_grid = Grids()

results = Initialize();

# First iteration
v_next, pf_next = Bellman(results);

# value function plot
plot(x = k_grid, results.value_function[:,1,5,1])
plot!(x = k_grid, v_next[:,1,5,1])

# policy function plot
plot(x = k_grid,  pf_next[:,1,5,1])
plot!(x = k_grid,  k_grid)

# update
results.value_function = v_next;

# second iteration
v_next, pf_next = Bellman(results);

# value functino plot
plot(x = k_grid, results.value_function[:,1,5,1])
plot!(x = k_grid, v_next[:,1,5,1])

# policy function plot
plot(x = k_grid, pf_next[:,1,5,1])
plot!(x = k_grid, k_grid)

# update
results.value_function = v_next;

# third iteration
v_next, pf_next = Bellman(results);

# value function plot
plot(x = k_grid, results.value_function[:,1,5,1])
plot!(x = k_grid, v_next[:,1,5,1])

# policy function plot
plot(x = k_grid, pf_next[:,1,5,1])
plot!(x = k_grid, k_grid)

# update
results.value_function = v_next;

# fourth iteration
v_next, pf_next = Bellman(results);

# value function plot
plot(x = k_grid, results.value_function[:,1,5,1])
plot!(x = k_grid, v_next[:,1,5,1])

# policy function plot
plot(x = k_grid, pf_next[:,1,5,1])
plot!(x = k_grid, k_grid)

####################################################################################################################
############################# test Solve_Bellman() #################################################################
####################################################################################################################

include("03_solve_value_function.jl");

@unpack k_grid = Grids();

results = Initialize();

@elapsed results = Solve_Bellman(results; progress = true)

# value function plot
plot(x = k_grid, results.value_function[:,1,10,:])
plot!(x = k_grid, results.value_function[:,2,10,:])

# policy function plot
plot(x = k_grid, results.policy_function[:,1,:,1])
plot!(x = k_grid, k_grid)

####################################################################################################################
############################# test simulate_capital_path() #########################################################
####################################################################################################################

include("04_solve_model.jl");

@unpack k_grid = Grids();

results = Initialize();

@elapsed results = Solve_Bellman(results)

results = simulate_capital_path(results);

plot(results.K)

####################################################################################################################
############################# test simulate_capital_path() #########################################################
####################################################################################################################

include("04_solve_model.jl");

@elapsed results = Solve_model()