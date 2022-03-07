# Alex von Hafften
# ECON 810: Advanced Macro
# PS 3 - Part 2 - Code for running analysis
# Professor Carter Braxton

using Plots, Statistics, DataFrames

cd("/Users/vonhafften/Documents/UW Madison/problem_sets/econ_810a/ps3/")

include("part_2_model.jl")
include("part_2_simulation.jl")

P = Primitives()
G = Grids()


###################################################################################################
# solve model
###################################################################################################

R = Initialize(0.4)
Solve!(R; progress = true)

# simulate model
S = Initialize_simulation(8000)
Simulate_model!(S, R; progress = true)

###################################################################################################
# figures
###################################################################################################

wages = f.(S.h) .* S.ω .* S.e

# histogram of assets
histogram(reshape(S.b, S.N*P.T), bins = 20, legend = false)
title!("Distribution of Assets")
xlabel!("Assets")
ylabel!("Count")
savefig("figure_3a.png")

# histogram of wages
histogram(reshape(wages, S.N*P.T), bins = 20, legend = false)
title!("Distribution of Wages")
xlabel!("Wages")
ylabel!("Count")
savefig("figure_3b.png")

# unemployment rate
sum(reshape(1 .- S.e, S.N*P.T))/(S.N*P.T)

# average assets over lifecycle
plot(mean(S.b, dims = 1)', legend = false)
title!("Average Assets over Lifecycle")
ylabel!("Assets")
xlabel!("Model Age")
savefig("figure_3d_1.png")

# average wages over lifecycle
plot(mean(wages, dims = 1)', legend = false)
title!("Average Wages over Lifecycle")
ylabel!("Wages")
xlabel!("Model Age")
savefig("figure_3d_2.png")

# average earnings gain while employed
df = DataFrame(wages = reshape(wages[:,2:end], S.N*(P.T-1)), 
               lag_wages = reshape(wages[:,1:(end-1)], S.N*(P.T-1)), 
               employed = reshape(S.e[:,2:end], S.N*(P.T-1)),
               lag_employed = reshape(S.e[:,1:(end-1)], S.N*(P.T-1)))

df = filter(:employed => ==(1), df)
df = filter(:lag_employed => ==(1), df)
df = filter(:lag_wages => >(0.0), df)

mean((df.wages .- df.lag_wages)./df.lag_wages)

# Earnings around job loss
earnings_jl = earnings_around_jl(S, R)

plot(-4:1:8, mean(earnings_jl; dims = 1)', legend = false)
xlabel!("Time around Job Loss")
ylabel!("Earnings")
title!("Average Earnings around Job Loss")
savefig("figure_3f.png")

# Consumption around job loss
consumption_jl = consumption_around_jl(S, R)

plot(-4:1:8, mean(consumption_jl; dims = 1)', legend = false)
xlabel!("Time around Job Loss")
ylabel!("Consumption")
title!("Average Consumption around Job Loss")
savefig("figure_3g.png")

# solve model with higher unemployment benefit
R_2 = Initialize(0.4 * 1.1)
Solve!(R_2; progress = true)
S_2 = Initialize_simulation(8000)
Simulate_model!(S_2, R_2; progress = true)

# Earnings around job loss
earnings_jl_2 = earnings_around_jl(S_2, R_2)
plot(-4:1:8, mean(earnings_jl_2; dims = 1)', label = "z = 0.44", legend = :bottomright)
plot!(-4:1:8, mean(earnings_jl; dims = 1)', label = "z = 0.4")
xlabel!("Time around Job Loss")
ylabel!("Earnings")
title!("Average Earnings around Job Loss")
savefig("figure_3h_1.png")

# Consumption around job loss
consumption_jl_2 = consumption_around_jl(S_2, R_2)
plot(-4:1:8, mean(consumption_jl_2; dims = 1)', label = "z = 0.44", legend = :bottomright)
plot!(-4:1:8, mean(consumption_jl; dims = 1)', label = "z = 0.4")
xlabel!("Time around Job Loss")
ylabel!("Consumption")
title!("Average Consumption around Job Loss")
savefig("figure_3h_2.png")

# unemployment rate
sum(reshape(1 .- S_2.e, S_2.N*P.T))/(S_2.N*P.T)