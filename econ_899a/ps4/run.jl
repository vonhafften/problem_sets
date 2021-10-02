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
k_guess_0 = 3.513795391470152
k_guess_1 = 4.598758704122076

# Initial guesses for initial and terminal steady state labor
l_guess_0 = 0.3468221868691036
l_guess_1 = 0.36525017624196743

Solve_transition(θ_0, θ_1, k_guess_0, k_guess_1, l_guess_0, l_guess_1;
                 progress = true)
