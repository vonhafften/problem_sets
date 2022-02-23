# Alex von Hafften 
# FIN 970: Asset Pricing
# HW 1 Problem 1
# Professor Ivan Shaliastovich

# This code tests the gmm estimation of a linear regression in ./gmm.jl.

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/fin_970/ps1")

using Distributions

include("gmm.jl")

# simulate data
n = 1000

X_1 = rand(Normal(0, 1), n)
X_2 = rand(Normal(0, 1), n)
X_3 = rand(Normal(0, 1), n)
epsilon = rand(Normal(0, 1), n)

Y = 1 .+ 2 .* X_1 .+ 3 .* X_2 .+ 4 .* X_3 + epsilon
X = hcat(ones(n), X_1, X_2, X_3)

# estimate gmm
gmm_estimate = gmm(Y, X, 3)