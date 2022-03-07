# Alex von Hafften
# ECON 810: Advanced Macro
# PS 4 - Part 2 - Code for running analysis
# Professor Carter Braxton

using Plots, Statistics, DataFrames

cd("/Users/vonhafften/Documents/UW Madison/problem_sets/econ_810a/ps4/")

include("part_2_model.jl")
include("part_2_simulation.jl")

P = Primitives()
G = Grids()


###################################################################################################
# solve model
###################################################################################################

R = Initialize()

Solve_terminal_period!(R)