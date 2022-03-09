# Alex von Hafften
# ECON 810A - Project
# Bewley with long-term illiquid Assets

using Parameters, Plots, Statistics

cd("/Users/vonhafften/Documents/UW Madison/problem_sets/econ_810a/project/")

include("model.jl")

P = Primitives()

# baseline
R_1 = Initialize(0.99, 0.2)
Solve!(R_1)

which_t = 28
plot(R_1.pf_m_l[which_t, 1, :, :]')

include("simulation.jl")
S_1 = Initialize_Simulation(10000)
Simulate!(S_1, R_1)

# policy experiment #1
R_2 = Initialize(0.95, 0.2)
Solve!(R_2)
S_2 = Initialize_Simulation(100)
Simulate!(S_2, R_2)

plot(mean(S_1.b; dims =1)', label = "β_I = 0.99")
plot!(mean(S_2.b; dims =1)', label = "β_I = 0.95")
title!("Average Bond Holdings over Lifecycle")
xlabel!("Model Age")
ylabel!("Bonds")
savefig("bond.png")

plot(mean(S_1.a; dims =1)', label = "β_I = 0.99")
plot!(mean(S_2.a; dims =1)', label = "β_I = 0.95")
title!("Average Cash Holdings over Lifecycle")
xlabel!("Model Age")
ylabel!("Cash")
savefig("cash.png")

plot(mean(S_1.b; dims =1)'./mean(S_1.a + S_1.b; dims =1)' .* 100, label = "β_I = 0.99")
plot!(mean(S_2.b; dims =1)'./mean(S_2.a + S_2.b; dims =1)' .* 100, label = "β_I = 0.95", legend = :bottomright)
title!("Average Portfolio Breakdown over Lifecycle")
xlabel!("Model Age")
ylabel!("Bond / (Cash + Bond) (%)")
savefig("portfolio.png")

plot(S_1.b[1:10, :]', legend = false)
title!("Bond Holdings for Select Simulations")
xlabel!("Model Age")
ylabel!("Bonds")
savefig("problem.png")