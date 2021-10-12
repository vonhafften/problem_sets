# Krusell and Smith (1998)
# Alex von Hafften
# October 13, 2021

# ECON 899A Computational Economics
# Problem set 5

# This file contains test functions.

using Plots

include("model.jl")

####################################################################################################################
####################################################################################################################
############################# test simulate_Z() ####################################################################
####################################################################################################################
####################################################################################################################

compute_Π_star(get_Π_z())

n_trials = 10000
recession_share = zeros(n_trials)

for i = 1:n_trials
    test = simulate_Z()
    recession_share[i] = count(test .== 0.99)/length(test)
end

histogram(recession_share) # should be centered at the stationary distribution