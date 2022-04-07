# Alex von Hafften
# Problem set 1
# ECON 717B: Applied Econometrics
# Professor Matt Wiswall
# April 7, 2022

# This file calls functions from ./analysis.jl

cd("/Users/vonhafften/Documents/UW Madison/problem_sets/econ_717b/ps1/")
include("analysis.jl")

using Plots

################################################################################################
# part a - why we can normalize π_1 and π_2
################################################################################################

# see write-up

################################################################################################
# part b - program to compute the model
################################################################################################

# see run.jl

################################################################################################
# part c - choose vector of parameters such that 60 percent of simulations choose occupation 1
################################################################################################

# true parameters
θ_0 = parameters(1.0, 1.0, 1.1, 1.0, 0.3, 0.5, 0.25)

# simulate true data
data_0 = simulate_data(θ_0;seed = 345)

# fraction that choose occupation is about 60 percent
sum(data_0.D)/data_0.N

################################################################################################
# part d - write mle estimator for θ
################################################################################################

# see write-up

################################################################################################
# part e - write code to estimate θ using mle
################################################################################################

(μ_1_hat, ρ_hat) = mle(θ_0, data_0)

################################################################################################
# part f - plot identification figure
################################################################################################

id_plot_μ_1(μ_1_hat, ρ_hat, θ_0, data_0)
id_plot_ρ(μ_1_hat, ρ_hat, θ_0, data_0)