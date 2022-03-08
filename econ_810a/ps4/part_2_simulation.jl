# Alex von Hafften
# ECON 810: Advanced Macro
# PS 4 - Part 2 - Code for simulating the model
# Professor Carter Braxton

using Parameters

cd("/Users/vonhafften/Documents/UW Madison/problem_sets/econ_810a/ps4/")

include("part_2_model.jl")

mutable struct Simulation_Results
    N::Int64                # number of simulations
    h::Array{Float64}       # human capital
    k::Array{Float64}       # asset holdings
    s::Array{Float64}       # time investment in human capital
    z::Array{Float64}       # log normal human capital shock
end

function Initialize_simulation(N::Int64)
    P = Primitives()
    G = Grids()

    h = zeros(N, P.T)
    k = zeros(N, P.T)
    s = zeros(N, P.T)
    z = zeros(N, P.T)

    Simulation_Results(N, h, k, s, z)
end

function Simulate_model!(S::Simulation_Results, R::Results)
    P = Primitives()
    G = Grids()

    # pull initial human capital level
    S.h[:, 1] = max.(min.(rand(R.distribution_h_0, S.N), G.max_h), G.min_h)

    for i = 1:S.N, t = 1:(P.T-1)

        # get asset level and index
        k = S.k[i,t]
        i_k = argmin(abs.(k .- G.grid_k))

        # get human capital level and index
        h = S.h[i,t]
        i_h = argmin(abs.(h .- G.grid_h))

        # makes time choice
        s = R.pf_s[t, i_h, i_k]
        S.s[i,t] = s

        # makes saving choice
        k_p = R.pf_k[t, i_h, i_k]
        S.k[i,t+1] = k_p

        # draw z shock
        z = rand(R.distribution_z)
        S.z[i, t+1] = z

        # update human capital
        S.h[i, t+1] = min(max(exp(z) * H(h, s, P.α), G.min_h), G.max_h)
    end
end

