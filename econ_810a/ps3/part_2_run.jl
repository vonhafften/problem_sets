# Alex von Hafften
# ECON 810: Advanced Macro
# PS 3 - Part 2 - Code for running analysis
# Professor Carter Braxton

using Plots

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_810a/ps3/")

include("part_2_model.jl")
include("part_2_simulation.jl")

P = Primitives()
G = Grids()

# solve model
R = Initialize()
Solve!(R; progress = true)

# simulate model
S = Initialize_simulation(100)
Simulate_model!(S, R; progress = true)

# figures

# histogram of assets
histogram(reshape(S.b, S.N*P.T))

# histogram of wages
histogram(reshape(f.(S.h) .* S.ω, S.N*P.T))

# unemployment rate
sum(reshape(1 .- S.e, S.N*P.T))/(S.N*P.T)