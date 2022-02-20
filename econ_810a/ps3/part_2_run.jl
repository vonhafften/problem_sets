# Alex von Hafften
# ECON 810: Advanced Macro
# PS 3 - Part 2 - Code for running analysis
# Professor Carter Braxton

using Plots

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_810a/ps3/")

include("part_2_model.jl")

# solve model
R = Initialize()
Solve!(R; progress = true)

# simulate model
