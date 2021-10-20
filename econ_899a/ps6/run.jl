####################################################################################################################
# ECON 899A Computational Economics
# Problem Set 6

# Hopenhayn and Rogerson (1993)
# Alex von Hafften
# October 20, 2021

# This file runs the model.
####################################################################################################################

using Plots, DataFrames

include("model.jl");

####################################################################################################################
########################################## Standard ################################################################
####################################################################################################################

@unpack ν = Primitives()

results = Initialize(1.0);
Solve_price(results);
Solve_M(results);

println("Print level: ", results.p)
println("Mass of incumbants: ", sum((1 .- results.x) .* results.μ))
println("Mass of entrants: ", results.M)
println("Mass of exits: ", sum(results.x .* results.μ))
println("Aggregate labor: ", results.L_d)
println("Labor of Incumbants: ", sum(results.N_d .* results.μ))
println("Labor of Entrants: ", results.M * sum(results.N_d .* ν))
println("Fraction of Labor in Entrants: ", results.M * sum(results.N_d .* ν)/ results.L_d)