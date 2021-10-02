# Conesa and Krueger (1999)
# Alex von Hafften
# October 6, 2021

# ECON 899A Computational Economics
# Problem Set 4

# This file calls transition.jl and runs the computations

using Plots

include("transition.jl");

################################################################################
######################## Exercise 1.1 ##########################################
################################################################################

# Initial and terminal steady states
θ_0 = 0.11
θ_1 = 0.0

# Initial guesses for initial and terminal steady state capital
k_0_0 = 3.518053688978396
k_0_1 = 4.596865227600611

# Initial guesses for initial and terminal steady state labor
l_0_0 = 0.34699359276158837
l_0_1 = 0.3651178450225594

Solve_transition(θ_0, θ_1, k_0_0, k_0_1, l_0_0, l_0_1)
