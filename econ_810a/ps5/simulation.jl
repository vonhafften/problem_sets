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

    h = zeros(N, 9)
    b = zeros(N, 9)
    i = zeros(N, 4)
    τ = zeros(N)

    for i = 1:N
        h[i, :] = simulate(P.tauchen_h, 9) # simulate human capital for whole life
    end

    b[:, 1] .= 100000.0

    Simulation_Results(N, h, b, i, τ)
end

function Simulate_generation!(R::Results, S::Simulation_Results)
    P = Primitives()

    for i = 1:S.N
        # t = 4
        i_b = argmin(abs.(S.b[i, 4-3] .- P.grid_b))
        i_h = argmin(abs.(S.h[i, 4-3] .- P.grid_h))
        S.b[i, 5-3] = R.pf_b_4[i_b, i_h]

        # t = 5
        h_c = sample(P.grid_h, Weights(P.Π_h_0), 1)[1] # draw child's starting human capital

        i_b = argmin(abs.(S.b[i, 5-3] .- P.grid_b))
        i_h = argmin(abs.(S.h[i, 5-3] .- P.grid_h))
        i_h_c = argmin(abs.(P.grid_h .- h_c))

        S.i[i, 6-5] = R.pf_i_5[i_b, i_h, i_h_c]
        S.b[i, 6-3] = R.pf_b_5[i_b, i_h, i_h_c]

        # t = 6
        h_c = log((1-P.ω_c)*exp(h_c) + P.γ*P.ω_c * S.i[i, 6-5])

        i_b = argmin(abs.(S.b[i, 6-3] .- P.grid_b))
        i_h = argmin(abs.(S.h[i, 6-3] .- P.grid_h))
        i_h_c = argmin(abs.(P.grid_h .- h_c))

        S.i[i, 7-5] = R.pf_i_6[i_b, i_h, i_h_c]
        S.b[i, 7-3] = R.pf_b_6[i_b, i_h, i_h_c]

        # t = 7
        h_c = log((1-P.ω_c)*exp(h_c) + P.γ*P.ω_c * S.i[i, 7-5])

        i_b = argmin(abs.(S.b[i, 7-3] .- P.grid_b))
        i_h = argmin(abs.(S.h[i, 7-3] .- P.grid_h))
        i_h_c = argmin(abs.(P.grid_h .- h_c))

        S.i[i, 8-5] = R.pf_i_7[i_b, i_h, i_h_c]
        S.b[i, 8-3] = R.pf_b_7[i_b, i_h, i_h_c]

        # t = 8
        h_c = log((1-P.ω_c)*exp(h_c) + P.γ*P.ω_c * S.i[i, 8-5])

        i_b = argmin(abs.(S.b[i, 8-3] .- P.grid_b))
        i_h = argmin(abs.(S.h[i, 8-3] .- P.grid_h))
        i_h_c = argmin(abs.(P.grid_h .- h_c))

        S.i[i, 9-5] = R.pf_i_8[i_b, i_h, i_h_c]
        S.b[i, 9-3] = R.pf_b_8[i_b, i_h, i_h_c]

        # t = 9
        h_c = log((1-P.ω_c)*exp(h_c) + P.γ*P.ω_c * S.i[i, 9-5])

        i_b = argmin(abs.(S.b[i, 9-3] .- P.grid_b))
        i_h = argmin(abs.(S.h[i, 9-3] .- P.grid_h))
        i_h_c = argmin(abs.(P.grid_h .- h_c))

        S.τ[i] = R.pf_τ_9[i_b, i_h, i_h_c]
        S.b[i, 10-3] = R.pf_b_9[i_b, i_h, i_h_c]

        # t = 10
        i_b = argmin(abs.(S.b[i, 10-3] .- P.grid_b))
        i_h = argmin(abs.(S.h[i, 10-3] .- P.grid_h))

        S.b[i, 11-3] = R.pf_b_10[i_b, i_h]

        # t = 11
        i_b = argmin(abs.(S.b[i, 11-3] .- P.grid_b))
        i_h = argmin(abs.(S.h[i, 11-3] .- P.grid_h))

        S.b[i, 12-3] = R.pf_b_11[i_b, i_h]
    end
end

function Simulate!(R::Results, S::Simulation_Results)
    P = Primitives()
    
    i = 1
    max_iter = 1000

    while true
        if i > max_iter
            break
        end
        
        println("Generation #", i)

        # randomly select simulation to be parent
        for i=1:S.N
            j = sample(1:S.N)

            j_b = argmin(abs.(S.b[j, 9-3] .- P.grid_b))
            j_h = argmin(abs.(S.h[j, 9-3] .- P.grid_h))
            i_h_c = argmin(abs.(S.h[i, 4-3] .- P.grid_h))
            
            S.b[i, 1] = R.pf_τ_9[j_b, j_h, i_h_c]
        end

        T_star!(R, S)

        i += 1
    end
end