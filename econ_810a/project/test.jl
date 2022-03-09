# Alex von Hafften
# ECON 810A - Project
# Bewley with long-term illiquid Assets

using Parameters, Plots

cd("/Users/vonhafften/Documents/UW Madison/problem_sets/econ_810a/project/")

include("model.jl")

P = Primitives()
R = Initialize(.99, 0.2)
Solve_terminal_period!(R)

R.vf_m[P.T, 1, :, :]
R.vf_m[P.T, 2, :, :]
R.vf_u[P.T, 1, :, :]
R.vf_u[P.T, 2, :, :]

plot(R.vf_m[P.T, 1, [1, P.N_a], :]', label=["low cash" "high cash"])
title!("Employed, matched, t = 30")
xlabel!("LT bond holdings")

plot(R.vf_m[P.T, 2, [1, P.N_a], :]', label=["low cash" "high cash"])
title!("Unemployed, matched, t = 30")
xlabel!("LT bond holdings")

plot(R.vf_m[P.T, 1, :, [1, P.N_b]], label=["low LT bonds" "high LT bonds"])
title!("Employed, matched, t = 30")
xlabel!("cash holdings")

plot(R.vf_m[P.T, 2, :, [1, P.N_b]], label=["low LT bonds" "high LT bonds"])
title!("Unemployed, matched, t = 30")
xlabel!("cash holdings")

plot(R.vf_u[P.T, 1, [1, P.N_a], :]')
xlabel!("Cash holdings")
plot(R.vf_u[P.T, 2, [1, P.N_a], :]')
xlabel!("Cash holdings")

Solve_nonterminal_period!(R)

R.vf_u[P.T-1, 1, :, :]
R.vf_u[P.T-1, 2, :, :]

plot(R.vf_u[P.T-1, 1, [1, 10, P.N_a-1], :]')
plot(R.vf_u[P.T-1, 2, [1, 10, P.N_a-1], :]')
plot(R.vf_u[P.T-1, 1, :, [1, 10, P.N_b-1]])
plot(R.vf_u[P.T-1, 2, :, [1, 10, P.N_b-1]])


plot(R.pf_u_b_tilde[P.T-1, 2, :, [1, 10, P.N_b]])
plot(R.pf_u_a[P.T-1, 2, :, [1, 10, P.N_b]])
plot(R.pf_u_c[P.T-1, 2, :, [1, 10, P.N_b]])

plot(R.vf_m[P.T-1, 1, [1, 10, P.N_a-1], :]')
plot(R.vf_m[P.T-1, 2, [1, 10, P.N_a-1], :]')
plot(R.vf_m[P.T-1, 1, :, [1, 10, P.N_b-1]])
plot(R.vf_m[P.T-1, 2, :, [1, 10, P.N_b-1]])

plot(R.pf_m_b_tilde[P.T-1, 2, :, [1, 10, P.N_b]])
plot(R.pf_m_a[P.T-1, 2, :, [1, 10, P.N_b]])
plot(R.pf_m_c[P.T-1, 2, :, [1, 10, P.N_b]])

which_t = 10

R.vf_u[which_t, 1, :, :]
R.vf_u[which_t, 2, :, :]

plot(R.vf_u[which_t, 1, [1, 10, P.N_a-1], :]')
plot(R.vf_u[which_t, 2, [1, 10, P.N_a-1], :]')
plot(R.vf_u[which_t, 1, :, [1, 10, P.N_b-1]])
plot(R.vf_u[which_t, 2, :, [1, 10, P.N_b-1]])

plot(R.pf_u_b_tilde[which_t, 1, :, [1, 10, P.N_b]])
plot(R.pf_u_b_tilde[which_t, 1, [1, 10, P.N_b], :]')
plot(R.pf_u_a[which_t, 1, :, [1, 10, P.N_b]])
plot(R.pf_u_a[which_t, 1, [1, 10, P.N_b], :]')
plot(R.pf_u_c[which_t, 1, :, [1, 10, P.N_b]])
plot(R.pf_u_c[which_t, 1, [1, 10, P.N_b], :]')

plot(R.pf_u_b_tilde[which_t, 2, :, [1, 10, P.N_b]])
plot(R.pf_u_a[which_t, 2, :, [1, 10, P.N_b]])
plot(R.pf_u_c[which_t, 2, :, [1, 10, P.N_b]])

plot(R.vf_m[which_t, 1, [1, 10, P.N_a-1], :]')
plot(R.vf_m[which_t, 2, [1, 10, P.N_a-1], :]')
plot(R.vf_m[which_t, 1, :, [1, 10, P.N_b-1]])
plot(R.vf_m[which_t, 2, :, [1, 10, P.N_b-1]])

plot(R.pf_m_b_tilde[which_t, 2, :, [1, 10, P.N_b]])
plot(R.pf_m_a[which_t, 2, :, [1, 10, P.N_b]])
plot(R.pf_m_c[which_t, 2, :, [1, 10, P.N_b]])
plot(R.pf_m_l[which_t, 2, :, [1, 10, P.N_b]])

plot(R.pf_m_l[which_t, 1, [12, 13, 14, 15], :]')

############


using Plots

R = Initialize()
Solve_HH_problem!(R; progress = true)

# plots 

@unpack grid_a, δ = Primitives()

plot(grid_a, R.vf, label = ["Employed" "Unemployed"], legend = :bottomright);
title!("Value Function");
xlabel!("Assets")
savefig("vf.png")

plot(grid_a, R.pf_c, label = ["Employed" "Unemployed"], legend = :bottomright);
title!("Consumption Policy Function");
xlabel!("Assets")
savefig("consumption_pf.png")

plot(grid_a, R.pf_a, label =["Employed" "Unemployed"], legend = :bottomright);
plot!(grid_a, grid_a, label = "45 degree line");
title!("Asset Policy Function");
xlabel!("Assets")
savefig("asset_pf.png")

plot(grid_a, R.pf_l, label =["Employed" "Unemployed"], legend = :bottomright);
title!("Liquiation Policy Function");
xlabel!("Assets")
savefig("liquidation_pf.png")

plot(grid_a, (1 .-R.pf_l).*[grid_a grid_a].*(1-δ), label =["Employed" "Unemployed"], legend = :bottomright);
plot!(grid_a, grid_a*(1-δ), label = "Total Long-Term Assets");
title!("Long-Term Assets Function");
xlabel!("Assets")
savefig("liquidated_asset_pf.png")


plot(grid_a[:,1], R.pf_c[:,1], legend = :bottomright, label = "consumption")
plot!(grid_a[:,1], R.pf_c[:,1] .+ R.pf_a[:,1], legend = :bottomright, label = "consumption + savings")
plot!(grid_a, grid_a .+ 2, label = "Total")
plot!(grid_a, grid_a .+ 2, label = "45 degree line")
title!("Consumption Policy Function");
xlabel!("Assets")