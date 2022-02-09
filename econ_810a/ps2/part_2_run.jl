# Ljungqvist and Sargent (1998)
# Alex von Hafften
# January 8, 2022

# ECON 810A Advanced Macro Theory
# Problem Set 2 - Part 2

using Plots

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_810a/ps2/")

include("part_2_model.jl")

@unpack grid_h = Primitives()

R = Initialize()
vfi!(R)

plot(grid_h, R.S_pf[336:3:360,:]', label = transpose((336:3:360)));
title!("Search Effort Policy Function")
savefig("search_pf.png")
plot(grid_h, R.RW_pf[30:60:330,:]', label =transpose((30:60:330)))
savefig("reservation_wage_pf.png")
plot(grid_h, R.RW_pf[336:6:360,:]', label =transpose((336:6:360)))
savefig("reservation_wage_pf_end.png")

#################################

S = Initialize_simulation(10000)

simulate_model!(S, R)

plot(mean(S.employment .* S.h; dims = 1)', label = "Employed");
plot!(mean((1 .- S.employment) .* S.h; dims = 1)', label = "Unemployed");
title!("Mean Human Capital")
savefig("mean_human_capital.png")

compute_wage_growth(S)

jl_data = earnings_around_jl(S)
plot(-6:24, 0.07.-mean(jl_data; dims = 1)', legend = false)
savefig("jl.png")
