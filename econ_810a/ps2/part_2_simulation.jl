# Ljungqvist and Sargent (1998)
# Alex von Hafften
# January 8, 2022

# ECON 810A Advanced Macro Theory
# Problem Set 2 - Part 2

# Load libraries
using Parameters, Distributions

include("part_2_model.jl")

# Results structure
mutable struct Simulation
    trials::Int64                   # Number of individuals simulated
    employment::Array{Int64}        # employment status
    s::Array{Float64}               # search effort
    rw::Array{Float64}              # reservation wage
    h::Array{Float64}               # human capital level
    w::Array{Float64}               # wage
end

function Initialize_simulation(trials)
    @unpack T = Primitives()

    employment = zeros(trials, T)
    s = zeros(trials, T)
    rw = zeros(trials, T)
    h = zeros(trials, T)
    w = zeros(trials, T)

    Simulation(trials, employment, s, rw, h, w)
end

function simulate_model!(S::Simulation, R::Results)
    @unpack grid_h, T, δ, ψ_e, ψ_u, Δ, max_h, min_h, μ_w, σ_w = Primitives()

    for i = 1:S.trials

        # individuals starts unemployment with uniform draw from hc grid
        S.employment[i, 1] = 0
        S.h[i,1] = sample(grid_h)
        S.w[i,1] = 0

        for t = 2:T
            # employed last period
            if S.employment[i, t-1] == 1
                
                # hc increases
                if rand(Uniform(0, 1)) > ψ_e
                    S.h[i,t] = min(S.h[i,t-1] + Δ, max_h)
                else # hc stay same
                    S.h[i,t] = S.h[i,t]
                end

                # stay employed
                if rand(Uniform(0, 1)) > δ
                    S.employment[i, t] = 1
                    S.w[i,t] = S.w[i,t-1]
                else # become unemployed
                    S.employment[i, t] = 0
                    S.w[i,t] = 0
                end

            else # if S.employment[i, t-1] == 0

                # hc decreases
                if rand(Uniform(0, 1)) > ψ_u
                    S.h[i,t] = max(S.h[i,t-1] - Δ, min_h)
                else # hc stay same
                    S.h[i,t] = S.h[i,t]
                end

                i_h = argmin(abs.(S.h[i,t] .- grid_h))

                S.s[i,t] = R.S_pf[t, i_h]
                S.rw[i,t] = R.RW_pf[t, i_h]

                # gets offer
                if rand(Uniform(0, 1)) < π(S.s[i,t])
                    w_offer = rand(Normal(μ_w, σ_w))

                    # accept offer
                    if w_offer > S.rw[i,t]
                        S.w[i,t] = w_offer
                        S.employment[i,t] = 1
                    else # if w_offer < rw[i,t] # reject offer
                        S.employment[i,t] = 0
                    end
                else
                    S.employment[i,t] = 0
                end
            end
        end
    end
end

function compute_wage_growth(S::Simulation)
    @unpack T = Primitives()

    wage = S.h .* S.w
    lag_wage = hcat(zeros(S.trials), wage[:, 1:(T-1)])
    wage_growth = reshape((wage .- lag_wage)./lag_wage, S.trials*T)

    # drops negative, infinite, and nans
    wage_growth = wage_growth[wage_growth .> 0]
    wage_growth = wage_growth[wage_growth .< 1]
    wage_growth = wage_growth[isfinite.(wage_growth)]

    return mean(wage_growth)

end

function earnings_around_jl(S::Simulation)
    @unpack T = Primitives()

    j = 1
    max_j = 1000000
    results = zeros(max_j, 31)

    for t = 25:(T-24), i = 1:S.trials
        if S.employment[i, t-1] == 1 & S.employment[i,t] == 0 & S.employment[i,t+1] == 0
            indexes = (t-6):(t+24)
            results[j, :] = S.h[i, indexes] .* S.w[i, indexes]
            j += 1
        end

        if j > max_j
            break
        end
    end

    results
end