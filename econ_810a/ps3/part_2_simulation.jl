# Alex von Hafften
# ECON 810: Advanced Macro
# PS 3 - Part 2 - Code for simulating the model
# Professor Carter Braxton

using Parameters

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_810a/ps3/")

include("part_2_model.jl")

mutable struct Simulation_Results
    N::Int64                # number of simulations
    e::Array{Int64}         # employment status
    h::Array{Float64}       # human capital
    b::Array{Float64}       # asset holdings
    ω::Array{Float64}       # wage piece rates
    θ::Array{Float64}       # market tightnesses
end

function Initialize_simulation(N::Int64)
    P = Primitives()
    G = Grids()

    e = fill(0, N, P.T)
    h = fill(G.min_h, N, P.T)
    b = zeros(N, P.T)
    ω = zeros(N, P.T)
    θ = zeros(N, P.T)

    Simulation_Results(N, e, h, b, ω, θ)
end

function Simulate_model!(S::Simulation_Results, R::Results; progress::Bool = false)
    P = Primitives()
    G = Grids()

    for i = 1:S.N, t = 1:(P.T-1)
        if S.e[i, t] == 0 # if worker is unemployed

            # get asset level and index
            b = S.b[i,t]
            i_b = argmin(abs.(b .- G.grid_b))

            # get human capital level and index
            h = S.h[i,t]
            i_h = argmin(abs.(h .- G.grid_h))

            # makes saving choice
            S.b[i,t+1] = R.u_b_pf[t, i_b, i_h]

            # makes piece rate choice
            S.ω[i,t+1]   = R.u_ω_pf[t, i_b, i_h]
            i_ω_p        = argmin(abs.(S.ω[i,t+1] .- G.grid_ω))

            # human capital shock
            if rand() < P.p_L
                i_h_p = max(i_h - 1, 1)
            else
                i_h_p = i_h
            end
            S.h[i,t+1] = G.grid_h[i_h_p]

            # market tightness
            S.θ[i, t+1] = R.θ[t+1, i_ω_p, i_h_p]

            # employment shock
            if rand() < p(S.θ[i, t+1], P.ζ)
                S.e[i, t+1] = 1
            else
                S.e[i, t+1] = 0
            end

        elseif S.e[i, t] == 1 # if worker is employed

            # get asset level and index
            b = S.b[i,t]
            i_b = argmin(abs.(b .- G.grid_b))

            # get human capital level and index
            h = S.h[i,t]
            i_h = argmin(abs.(h .- G.grid_h))

            # get piece rate level and index
            ω = S.ω[i,t]
            i_ω = argmin(abs.(ω .- G.grid_ω))

            # makes saving choice
            S.b[i,t+1] = R.w_b_pf[t, i_ω, i_b, i_h]

            # human capital shock
            if rand() < P.p_H
                i_h_p = min(i_h + 1, G.N_h)
            else
                i_h_p = i_h
            end
            S.h[i,t+1] = G.grid_h[i_h_p]

            # employment shock
            if rand() < P.δ
                S.e[i, t+1] = 0
            else
                S.e[i, t+1] = 1
            end

        else
            error("Error in simulation")
        end
    end
    S
end