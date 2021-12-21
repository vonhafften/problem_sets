# Computational Economics
# Professor JF Houde
# Problem set 1
# Alex von Hafften, Michael Nattinger, Sarah Bass, Xinxin Hu
# November 7, 2021

using Plots, Tables, DataFrames, CSV, StatFiles, Optim

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_899b/ps1/")

include("data.jl");
include("model.jl");

##########################################################################################
# Problem 1
β_1 = prepend!(zeros(K-1), -1.0);

p1_ll = log_likelihood(β_1, x, y)

p1_score = score(β_1, x, y)
CSV.write("tables/p1_score.csv", Tables.table(p1_score))

p1_hessian = hessian(β_1, x, y)
CSV.write("tables/p1_hessian.csv", Tables.table(p1_hessian))

##########################################################################################
# Exercise 2
p2_score = score_numerical(β_1, x, y)
CSV.write("tables/p2_score.csv", Tables.table(p2_score))

p2_hessian = hessian_numerical(β_1, x, y)
CSV.write("tables/p2_hessian.csv", Tables.table(p2_hessian))

##########################################################################################
# Exercise 3
results_newton = newton(log_likelihood, score, hessian, β_1, x, y)
CSV.write("tables/p3_beta.csv", Tables.table(results_newton[1]))

##########################################################################################
# Exercise 4

# Method without derivatives
results_simplex = optimize(β -> -log_likelihood(β, x, y), β_1, NelderMead())                 # fails if guess is β_1
results_simplex_2 = optimize(β -> -log_likelihood(β, x, y), results_newton[1], NelderMead()) # it does work if guess is newton results

# Method with gradient
results_lbfgs = optimize(β -> -log_likelihood(β, x, y), score!, β_1, LBFGS())

# (Built-in) with gradient and hessian
results_newton_2 = optimize(β -> -log_likelihood(β, x, y), score!, hessian!, β_1, Newton())

comparison_table  = Tables.table(hcat(results_simplex.minimizer, results_simplex_2.minimizer, results_lbfgs.minimizer, results_newton_2.minimizer))

CSV.write("tables/p4_beta_comparison.csv", comparison_table)
