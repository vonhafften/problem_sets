# Alex von Hafften
# ECON 810: Advanced Macro
# PS 3 - Part 2 - Code for testing the model code
# Professor Carter Braxton

using Plots

cd("/Users/vonhafften/Documents/UW Madison/problem_sets/econ_810a/ps3/")

include("part_2_model.jl")

P = Primitives()
G = Grids()

###################################################################################################
###################################################################################################

# test p and p_f and p_f_inv

θs = 0:0.01:10

job_finding_rates = p.(θs, P.ζ)
plot(θs, job_finding_rates)

hiring_rates = p_f.(θs, P.ζ)
plot(θs, hiring_rates)

# should equal θs
p_f_inv.(hiring_rates, P.ζ)

js = 0:0.1:10
θs = free_entry.(js, P.κ, P.ζ)
plot(js, θs)

###################################################################################################
###################################################################################################

# test Solve_terminal_period

R = Initialize()
Solve_terminal_period!(R; progress = true)

plot(G.grid_ω, R.j_vf[P.T, :, :], label = G.grid_h')
title!("Terminal Period: Firm value over hc levels")
ylabel!("firm value")
xlabel!("wage piecerate")

plot(G.grid_ω, R.θ[P.T, :, :], label = G.grid_h')
title!("Terminal Period: Market tightness over hc levels")
ylabel!("Market tightness")
xlabel!("wage piecerate")

plot(G.grid_b, R.u_vf[P.T, :, :], label = G.grid_h')
title!("Terminal Period: Unemployeds value over hc levels")
ylabel!("Unemployeds value")
xlabel!("asset holding")

plot(G.grid_b[2:G.N_b], R.w_vf[P.T, :, 2:G.N_b, 1]', label = G.grid_ω')
title!("Terminal Period: Employeds value over piecerates (min hc)")
ylabel!("Employeds value")
xlabel!("asset holding")

plot(G.grid_b[2:G.N_b], R.w_vf[P.T, :, 2:G.N_b, G.N_h]', label = G.grid_ω')
title!("Terminal Period: Employeds value over piecerates (max hc)")
ylabel!("Employeds value")
xlabel!("asset holding")


###################################################################################################
###################################################################################################

# test Solve_nonterminal_period

Solve_nonterminal_periods!(R; progress = true)

# Penultimate period

plot(G.grid_ω, R.j_vf[P.T-1, :, :], label = G.grid_h')
title!("Penultimate Period: Firm value over hc levels")
ylabel!("firm value")
xlabel!("wage piecerate")

plot(G.grid_ω, R.θ[P.T-1, :, :], label = G.grid_h')
title!("Penultimate Period: Market tightness over hc levels")
ylabel!("Market tightness")
xlabel!("wage piecerate")

plot(G.grid_b, R.u_vf[P.T-1, :, :], label = G.grid_h')
title!("Penultimate Period: Unemployeds value over hc levels")
ylabel!("Unemployeds value")
xlabel!("asset holding")

plot(G.grid_b, R.u_b_pf[P.T-1, :, :], label = G.grid_h')
title!("Penultimate Period: Unemployeds asset policy over hc levels")
ylabel!("asset holding tomorrow")
xlabel!("asset holding today")

plot(G.grid_b, R.u_ω_pf[P.T-1, :, :], label = G.grid_h')
title!("Penultimate Period: Unemployeds search policy over hc levels")
ylabel!("piece rate tomorrow")
xlabel!("asset holding")

plot(G.grid_b[2:G.N_b], R.w_vf[P.T-1, :, 2:G.N_b, 1]', label = G.grid_ω')
title!("Penultimate Period: Employeds value over piecerates (min hc)")
ylabel!("Employeds value")
xlabel!("asset holding")

plot(G.grid_b[2:G.N_b], R.w_vf[P.T-1, :, 2:G.N_b, G.N_h]', label = G.grid_ω')
title!("Penultimate Period: Employeds value over piecerates (max hc)")
ylabel!("Employeds value")
xlabel!("asset holding")

plot(G.grid_b, R.w_b_pf[P.T-1, :, :, 1]', label = G.grid_ω')
title!("Penultimate Period: Employeds asset policy over piecerates (max hc)")
ylabel!("asset holdings tomorrow")
xlabel!("asset holding")

plot(G.grid_b, R.w_b_pf[P.T-1, :, :, G.N_h]', label = G.grid_ω')
title!("Penultimate Period: Employeds asset policy over piecerates (max hc)")
ylabel!("asset holdings tomorrow")
xlabel!("asset holding")

# first period

plot(G.grid_ω, R.j_vf[1, :, :], label = G.grid_h')
title!("First Period: Firm value over hc levels")
ylabel!("firm value")
xlabel!("wage piecerate")

plot(G.grid_ω, R.θ[1, :, :], label = G.grid_h')
title!("Penultimate Period: Market tightness over hc levels")
ylabel!("Market tightness")
xlabel!("wage piecerate")

plot(G.grid_b, R.u_vf[1, :, :], label = G.grid_h')
title!("first Period: Unemployeds value over hc levels")
ylabel!("Unemployeds value")
xlabel!("asset holding")

plot(G.grid_b, R.u_b_pf[1, :, :], label = G.grid_h')
title!("first Period: Unemployeds asset policy over hc levels")
ylabel!("asset holding tomorrow")
xlabel!("asset holding today")

plot(G.grid_b, R.u_ω_pf[1, :, :], label = G.grid_h')
title!("first Period: Unemployeds search policy over hc levels")
ylabel!("piece rate tomorrow")
xlabel!("asset holding")

plot(G.grid_b[2:G.N_b], R.w_vf[1, :, 2:G.N_b, 1]', label = G.grid_ω')
title!("First Period: Employeds value over piecerates (min hc)")
ylabel!("Employeds value")
xlabel!("asset holding")

plot(G.grid_b[2:G.N_b], R.w_vf[1, :, 2:G.N_b, G.N_h]', label = G.grid_ω')
title!("First Period: Employeds value over piecerates (max hc)")
ylabel!("Employeds value")
xlabel!("asset holding")

plot(G.grid_b, R.w_b_pf[1, :, :, 1]', label = G.grid_ω')
title!("First Period: Employeds asset policy over piecerates (max hc)")
ylabel!("asset holdings tomorrow")
xlabel!("asset holding")

plot(G.grid_b, R.w_b_pf[1, :, :, G.N_h]', label = G.grid_ω')
title!("First Period: Employeds asset policy over piecerates (max hc)")
ylabel!("asset holdings tomorrow")
xlabel!("asset holding")

###################################################################################################
###################################################################################################

# test simulation code

S = Initialize_simulation(1)
Simulate_model!(S, R; progress = true)

plot(S.e')
plot(S.h')
plot(S.b')
plot(S.ω')
plot(S.θ')