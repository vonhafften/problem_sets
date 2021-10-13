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

n_trials = 10000
recession_share = zeros(n_trials)

for i = 1:n_trials
    test = simulate_Z()
    recession_share[i] = count(test .== 0.99)/length(test)
end

compute_Π_star(get_Π_z())

histogram(recession_share) # should be centered at the stationary distribution

sum(recession_share)/length(recession_share)


####################################################################################################################
####################################################################################################################
############################# test simulate_E(Z) ###################################################################
####################################################################################################################
####################################################################################################################

n_trials = 100;

recession_employment = zeros(n_trials);
recession_unemployment = zeros(n_trials);
boom_employment = zeros(n_trials);
boom_unemployment = zeros(n_trials);

# takes a bit...
for i = 1:n_trials
    println(i)
    Z = simulate_Z()
    test = simulate_E(Z)

    test_recession = test[Z .== 0.99, :]
    test_boom= test[Z .== 1.01, :]

    recession_employment[i] = count(test_recession .== 1)/length(test)
    boom_employment[i] = count(test_boom .== 1)/length(test)
    recession_unemployment[i] = count(test_recession .== 0)/length(test)
    boom_unemployment[i] = count(test_boom .== 0)/length(test)
end

compute_Π_star(get_Π_ε())

histogram(recession_employment);
histogram(boom_employment);
histogram(recession_unemployment);
histogram(boom_unemployment);

sum(boom_employment)/length(boom_employment)
sum(recession_employment)/length(recession_employment)
sum(boom_unemployment)/length(boom_unemployment)
sum(recession_unemployment)/length(recession_unemployment)


####################################################################################################################
####################################################################################################################
############################# test Initialize() ######## ###########################################################
####################################################################################################################
####################################################################################################################

results = Initialize()
