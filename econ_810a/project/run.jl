# Alex von Hafften
# ECON 810A - Project
# Bewley with long-term illiquid Assets

using Parameters, Plots, Statistics

cd("/Users/vonhafften/Documents/UW Madison/problem_sets/econ_810a/project/")


include("model.jl")
include("simulation.jl")

P = Primitives()

# baseline
R_1 = Solve(0.99, 0.2)
S_1 = Initialize_Simulation(10000)
Simulate!(S_1, R_1)


plot(mean(S_1.b; dims =1)', label = "β_I = 0.99 and γ = 0.2")
title!("Average Bond Holdings over Lifecycle")
xlabel!("Model Age")
ylabel!("Bonds")
savefig("bond_bl.png")

plot(mean(S_1.a; dims =1)', label = "β_I = 0.99 and γ = 0.2")
title!("Average Cash Holdings over Lifecycle")
xlabel!("Model Age")
ylabel!("Cash")
savefig("cash_bl.png")

plot(mean(S_1.b; dims =1)'./mean(S_1.a + S_1.b; dims =1)' .* 100, label = "β_I = 0.99 and γ = 0.2", legend = :bottomright)
title!("Average Portfolio Breakdown over Lifecycle")
xlabel!("Model Age")
ylabel!("Bond / (Cash + Bond) (%)")
savefig("portfolio_bl.png")


# policy experiment #1
R_2 = Solve(0.95, 0.2)
S_2 = Initialize_Simulation(10000)
Simulate!(S_2, R_2)

# policy experiment #2
R_3 = Solve(0.99, 0.1)
S_3 = Initialize_Simulation(10000)
Simulate!(S_3, R_3)

# plots

plot(mean(S_1.b; dims =1)', label = "β_I = 0.99 and γ = 0.2")
plot!(mean(S_2.b; dims =1)', label = "β_I = 0.95 and γ = 0.2")
plot!(mean(S_3.b; dims =1)', label = "β_I = 0.99 and γ = 0.1")
title!("Average Bond Holdings over Lifecycle")
xlabel!("Model Age")
ylabel!("Bonds")
savefig("bond.png")

plot(mean(S_1.a; dims =1)', label = "β_I = 0.99 and γ = 0.2")
plot!(mean(S_2.a; dims =1)', label = "β_I = 0.95 and γ = 0.2")
plot!(mean(S_3.a; dims =1)', label = "β_I = 0.99 and γ = 0.1")
title!("Average Cash Holdings over Lifecycle")
xlabel!("Model Age")
ylabel!("Cash")
savefig("cash.png")

plot(mean(S_1.b; dims =1)'./mean(S_1.a + S_1.b; dims =1)' .* 100, label = "β_I = 0.99 and γ = 0.2", legend = :bottomright)
plot!(mean(S_2.b; dims =1)'./mean(S_2.a + S_2.b; dims =1)' .* 100, label = "β_I = 0.95 and γ = 0.2")
plot!(mean(S_3.b; dims =1)'./mean(S_3.a + S_3.b; dims =1)' .* 100, label = "β_I = 0.99 and γ = 0.1")
title!("Average Portfolio Breakdown over Lifecycle")
xlabel!("Model Age")
ylabel!("Bond / (Cash + Bond) (%)")
savefig("portfolio.png")
