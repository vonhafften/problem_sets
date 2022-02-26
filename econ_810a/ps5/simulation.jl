# Alex von Hafften
# ECON 810: Advanced Macro
# PS 5 - Code for simulating the model
# Professor Carter Braxton

cd("/Users/vonhafften/Documents/uw_madison/problem_sets/econ_810a/ps5/")

include("model.jl")

using QuantEcon, StatsBase

# results
mutable struct Simulation_Results
    N::Int64             # number of trials
    h::Matrix{Float64}   # human capital
    h_c::Matrix{Float64} # human capital of child
    b::Matrix{Float64}   # asset holdings
    i::Matrix{Float64}   # investment in child
    τ::Vector{Float64}   # transfer to child
end

function Initialize_simulation(N::Int64)
    P = Primitives()

    h   = zeros(N, 9)
    h_c = zeros(N, 5)
    b   = zeros(N, 9)
    i   = zeros(N, 4)
    τ   = zeros(N)

    # set initial simulation set
    h[:, 1] .= P.min_h
    b[:, 1] .= 100000.0

    Simulation_Results(N, h, h_c, b, i, τ)
end

function Simulate_generation!(R::Results, S::Simulation_Results)
    P = Primitives()

    for i = 1:S.N
        # t = 4
        i_b = argmin(abs.(S.b[i, 4-3] .- P.grid_b))
        i_h = argmin(abs.(S.h[i, 4-3] .- P.grid_h))

        S.b[i, 5-3] = R.pf_b_4[i_b, i_h]
        S.h[i, 5-3] = sample(P.grid_h, Weights(P.Π_h[i_h, :]), 1)[1]

        # t = 5
        S.h_c[i, 5-4] = sample(P.grid_h, Weights(P.Π_h_0), 1)[1] # draw child's starting human capital

        i_b = argmin(abs.(S.b[i, 5-3] .- P.grid_b))
        i_h = argmin(abs.(S.h[i, 5-3] .- P.grid_h))
        i_h_c = argmin(abs.(S.h_c[i, 5-4] .- P.grid_h))

        S.i[i, 6-5] = R.pf_i_5[i_b, i_h, i_h_c]
        S.b[i, 6-3] = R.pf_b_5[i_b, i_h, i_h_c]
        S.h[i, 6-3] = sample(P.grid_h, Weights(P.Π_h[i_h, :]), 1)[1]

        # t = 6
        S.h_c[i, 6-4] = log((1-P.ω_c)*exp(S.h_c[i, 5-4]) + P.γ*P.ω_c * S.i[i, 6-5])

        i_b = argmin(abs.(S.b[i, 6-3] .- P.grid_b))
        i_h = argmin(abs.(S.h[i, 6-3] .- P.grid_h))
        i_h_c = argmin(abs.(S.h_c[i, 6-4] .- P.grid_h))

        S.i[i, 7-5] = R.pf_i_6[i_b, i_h, i_h_c]
        S.b[i, 7-3] = R.pf_b_6[i_b, i_h, i_h_c]
        S.h[i, 7-3] = sample(P.grid_h, Weights(P.Π_h[i_h, :]), 1)[1]

        # t = 7
        S.h_c[i, 7-4] = log((1-P.ω_c)*exp(S.h_c[i, 6-4]) + P.γ*P.ω_c * S.i[i, 7-5])

        i_b = argmin(abs.(S.b[i, 7-3] .- P.grid_b))
        i_h = argmin(abs.(S.h[i, 7-3] .- P.grid_h))
        i_h_c = argmin(abs.(S.h_c[i, 7-4] .- P.grid_h))

        S.i[i, 8-5] = R.pf_i_7[i_b, i_h, i_h_c]
        S.b[i, 8-3] = R.pf_b_7[i_b, i_h, i_h_c]
        S.h[i, 8-3] = sample(P.grid_h, Weights(P.Π_h[i_h, :]), 1)[1]

        # t = 8
        S.h_c[i, 8-4] = log((1-P.ω_c)*exp(S.h_c[i, 7-4]) + P.γ*P.ω_c * S.i[i, 8-5])

        i_b = argmin(abs.(S.b[i, 8-3] .- P.grid_b))
        i_h = argmin(abs.(S.h[i, 8-3] .- P.grid_h))
        i_h_c = argmin(abs.(S.h_c[i, 8-4] .- P.grid_h))

        S.i[i, 9-5] = R.pf_i_8[i_b, i_h, i_h_c]
        S.b[i, 9-3] = R.pf_b_8[i_b, i_h, i_h_c]
        S.h[i, 9-3] = sample(P.grid_h, Weights(P.Π_h[i_h, :]), 1)[1]

        # t = 9
        S.h_c[i, 9-4] = log((1-P.ω_c)*exp(S.h_c[i, 8-4]) + P.γ*P.ω_c * S.i[i, 9-5])

        i_b = argmin(abs.(S.b[i, 9-3] .- P.grid_b))
        i_h = argmin(abs.(S.h[i, 9-3] .- P.grid_h))
        i_h_c = argmin(abs.(S.h_c[i, 9-4] .- P.grid_h))

        S.τ[i] = R.pf_τ_9[i_b, i_h, i_h_c]
        S.b[i, 10-3] = R.pf_b_9[i_b, i_h, i_h_c]
        S.h[i, 10-3] = sample(P.grid_h, Weights(P.Π_h[i_h, :]), 1)[1]

        # t = 10
        i_b = argmin(abs.(S.b[i, 10-3] .- P.grid_b))
        i_h = argmin(abs.(S.h[i, 10-3] .- P.grid_h))

        S.b[i, 11-3] = R.pf_b_10[i_b, i_h]
        S.h[i, 11-3] = sample(P.grid_h, Weights(P.Π_h[i_h, :]), 1)[1]

        # t = 11
        i_b = argmin(abs.(S.b[i, 11-3] .- P.grid_b))
        i_h = argmin(abs.(S.h[i, 11-3] .- P.grid_h))

        S.b[i, 12-3] = R.pf_b_11[i_b, i_h]
        S.h[i, 12-3] = sample(P.grid_h, Weights(P.Π_h[i_h, :]), 1)[1]
    end
end

function Simulate!(R::Results, S::Simulation_Results)
    P = Primitives()
    
    i = 1
    max_iter = 200

    while i < max_iter
        
        println("Generation #", i)

        Simulate_generation!(R, S)

        S.h[:, 1] .= S.h_c[:, end]
        S.b[:, 1] .= S.τ

        i += 1
    end
end