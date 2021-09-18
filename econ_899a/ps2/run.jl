# Huggett (1993)
# Alex von Hafften
# September 16, 2021

# ECON 899A Computational Economics
# Problem set 2

include("model.jl");

using Plots

@unpack a_grid, α, β = Primitives();

################################################################################
##################################### Part II ##################################
################################################################################

results = Initialize();

Solve_model(results)

# value function at equilibrium price
plot(a_grid, results.value_function, labels = ["Employed" "Unemployed"])

# policy function at equilibrium price
plot(a_grid, [results.policy_function a_grid],
     labels = ["Employed" "Unemployed" "45° Line"],
     legend=:bottomright)
savefig("policy_function.png")

# a_hat is where the employed policy function crosses 45 degree line
# Dean's a_hat is around 1.0381
a_hat_e = a_grid[argmin(abs.(results.policy_function[:, 1] - a_grid))]
a_hat_u = a_grid[argmin(abs.(results.policy_function[:, 2] - a_grid))]

# Bond holding distribution
plot(a_grid, results.μ, labels = ["Employed" "Unemployed"])

# Calculate wealth distribution
w = calculate_wealth_distribution(results)
plot(a_grid, w, labels = ["Employed" "Unemployed"])
savefig("wealth_distribution.png")

# Where is the biggest wealth spike?
# Dean's w_spike is around 2.0381
w_spike = a_grid[argmax(w)]

# Calculate Lorenz curve
l = calculate_lorenz_curve(w);
plot([l[:,1] l[:,1]], [l[:,1] l[:,2]], labels = ["45° Line" "Lorenz Curve"],
     legend=:bottomright)
savefig("lorenz_curve.png")

gini_efficient = calculate_gini(l)

################################################################################
##################################### Part III #################################
################################################################################

# welfare in complete markets
# fb for first best
w_fb = calculate_w_fb()

# consumption equivalent
λ = calculate_λ(results, w_fb)
plot(a_grid, λ, labels = ["Employed" "Unemployed"])
savefig("lambda.png")

# welfare in incomplete markets
w_inc = sum(results.μ .* results.value_function)

# welfare gain
wg = sum(results.μ .* λ)

# fraction of population that prefers complete market
f = sum((λ .>= 0) .*results.μ)
