# Alex von Hafften
# ECON 810: Advanced Macro
# PS 4 - Part 2 - Code for running analysis
# Professor Carter Braxton

using Plots

cd("/Users/vonhafften/Documents/UW Madison/problem_sets/econ_810a/ps4/")

include("part_2_model.jl")
include("part_2_simulation.jl")

P = Primitives()
G = Grids()

# estimate and simulate baseline model

R_1 = Initialize(-0.029, 0.11, 2.0, 0.5)
Solve!(R_1)
S_1 = Initialize_simulation(100000)
Simulate_model!(S_1, R_1)

##################################################################################################################
# question 1 - earnings over lifecycle - mean, std dev, kurtosis
##################################################################################################################

earnings = S_1.h .* (1 .- S_1.s) .* P.R'

plot(mean(earnings; dims= 1)', legend = false)
title!("Mean Earnings over Lifecycle")
xlabel!("Model Age")
savefig("figure_1_mean.png")

plot(std(earnings; dims= 1)', legend = false)
title!("Standard Deviation of Earnings over Lifecycle")
xlabel!("Model Age")
savefig("figure_1_sd.png")

e_skewness = zeros(P.T)
for t in 1:P.T 
    e_skewness[t] = skewness(earnings[:, t])
end
plot(e_skewness, legend = false)
title!("Skewness of Earnings over Lifecycle")
xlabel!("Model Age")
savefig("figure_1_skewness.png")

e_kurtosis = zeros(P.T)
for t in 1:P.T 
    e_kurtosis[t] = kurtosis(earnings[:, t])
end
plot(e_kurtosis, legend = false)
title!("Kurtosis of Earnings over Lifecycle")
xlabel!("Model Age")
savefig("figure_1_kurtosis.png")

##################################################################################################################
# part b - policy functions over k and h
##################################################################################################################

which_t = 10

plot(G.grid_k, R_1.pf_s[which_t, [10, 50, 90], :]', labels = ["Low h" "Med. h" "High h"], legend = :bottomright)
xlabel!("Assets")
title!("Policy for Investing in Human Capital (t = 10)")
savefig("figure_2_10_pf_k.png")

plot(G.grid_h, R_1.pf_s[which_t, :, [10, 50, 90]], labels = ["Low k" "Med. k" "High k"])
xlabel!("Human Capital")
title!("Policy for Investing in Human Capital (t = 10)")
savefig("figure_2_10_pf_h.png")

which_t = 20

plot(G.grid_k, R_1.pf_s[which_t, [10, 50, 90], :]', labels = ["Low h" "Med. h" "High h"], legend = :bottomright)
xlabel!("Assets")
title!("Policy for Investing in Human Capital (t = 20)")
savefig("figure_2_20_pf_k.png")

plot(G.grid_h, R_1.pf_s[which_t, :, [10, 50, 90]], labels = ["Low k" "Med. k" "High k"])
xlabel!("Human Capital")
title!("Policy for Investing in Human Capital (t = 20)")
savefig("figure_2_20_pf_h.png")

##################################################################################################################
# part c - policy function with higher human capital shock variance
##################################################################################################################

R_2 = Initialize(-0.029, 0.22, 2.0, 0.5)
Solve!(R_2)

which_t = 10

plot(G.grid_k, R_2.pf_s[which_t, [10, 50, 90], :]', labels = ["Low h" "Med. h" "High h"], legend = :bottomright)
xlabel!("Assets")
title!("Policy for Investing in Human Capital (t = 10)")
savefig("figure_3_10_pf_k.png")

plot(G.grid_h, R_2.pf_s[which_t, :, [10, 50, 90]], labels = ["Low k" "Med. k" "High k"])
xlabel!("Human Capital")
title!("Policy for Investing in Human Capital (t = 10)")
savefig("figure_3_10_pf_h.png")

which_t = 20

plot(G.grid_k, R_2.pf_s[which_t, [10, 50, 90], :]', labels = ["Low h" "Med. h" "High h"], legend = :bottomright)
xlabel!("Assets")
title!("Policy for Investing in Human Capital (t = 20)")
savefig("figure_3_20_pf_k.png")

plot(G.grid_h, R_2.pf_s[which_t, :, [10, 50, 90]], labels = ["Low k" "Med. k" "High k"])
xlabel!("Human Capital")
title!("Policy for Investing in Human Capital (t = 20)")
savefig("figure_3_20_pf_h.png")

##################################################################################################################
# part d - increase sd of initial human capital distribution
##################################################################################################################

R_3 = Initialize(-0.029, 0.11, 2.0, 1.0)
Solve!(R_3)
S_3 = Initialize_simulation(100000)
Simulate_model!(S_3, R_3)

earnings_3 = S_3.h .* (1 .- S_3.s) .* P.R'

gksw_1 = mean(earnings; dims = 2 )
gksw_3 = mean(earnings_3; dims = 2)

histogram(gksw_1,alpha = .5, label = "Low σ_0")
histogram!(gksw_3,alpha = .5, label = "High σ_0", legend = :topleft)
title!("GKSW Lifetime Earnings")
xlabel!("Earnings")
ylabel!("Count")
savefig("figure_4_histogram")

plot(std(earnings; dims= 1)', label = "Low σ_0")
plot!(std(earnings_3; dims= 1)',  label = "High σ_0")
title!("Standard Deviation of Earnings over Lifecycle")
xlabel!("Model Age")
savefig("figure_4_sd.png")