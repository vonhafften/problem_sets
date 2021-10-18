####################################################################################################################
# Krusell and Smith (1998)
# Alex von Hafften
# October 13, 2021

# ECON 899A Computational Economics
# Problem set 5

# This file contains runs the model.
####################################################################################################################

include("04_solve_model.jl");

@elapsed results = Solve_model()

# print results for write-update

println(results.a0)
println(results.a1)
println(results.b0)
println(results.b1)
println(results.R2)