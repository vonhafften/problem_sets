# Alex von Hafften
# ECON 810: Advanced Macro
# PS 5 - Code for running the model
# Professor Carter Braxton

using Plots

cd("/Users/vonhafften/Documents/uw_madison/problem_sets/econ_810a/ps5/")

include("model.jl")
include("simulation.jl")

P = Primitives()
R = Initialize()

Solve!(R)

S = Initialize_simulation(10000)
Simulate!(R, S)

plot(4:12, mean(S.b, dims = 1)')
plot(4:12, var(S.b, dims = 1)')


