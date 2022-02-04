# Kaplan and Violante (2010)
# Alex von Hafften
# January 31, 2022

# ECON 810A Advanced Macro Theory
# Problem set 1 - Part 2

# This file contains the functions to estimate the model.
# ./part_2_run.jl calls this file and runs the model.

# Load libraries
using Parameters, QuantEcon, CSV, DataFrames

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_810a/ps1/")

# Primitive structure
@with_kw struct Primitives

    # life-cycle
    N::Int64                  = 35                               # lifespan
    κ::Array{Float64, 1}      = CSV.read("age_profile.csv", DataFrame).log_income # age profile of earnings from R code (in log income terms)

    # preference parameters
    γ::Float64                = 2.0                              # coefficient of relative risk aversion
    β::Float64                = 0.975                            # Discount rate

    # asset grid
    min_a::Float64            = 0.0                               # Lower bound
    max_a::Float64            = 2000000.0                         # Upper bound
    N_a::Int64                = 500                               # Number of points
    grid_srl_a                = range(min_a, max_a; length = N_a) # Step range
    grid_a::Array{Float64, 1} = collect(grid_srl_a)               # Grid array

    # transitory income shocks
    N_ε                       = 5                                 # Number states
    ρ_ε                       = 0.0                               # Persistence
    σ_ε                       = sqrt(0.05627207)                  # std dev - From R code
    mc_ε                      = tauchen(N_ε, ρ_ε, σ_ε)            # markov chain object
    grid_ε                    = mc_ε.state_values                 # states
    Π_ε                       = mc_ε.p                            # transition matrix

    # persistent income shocks
    N_ζ                       = 5                                 # Number of states
    ρ_ζ                       = 0.97                              # Persistence 
    σ_ζ                       = sqrt(0.06727484)                  # std dev - From R code
    mc_ζ                      = tauchen(N_ζ, ρ_ζ, σ_ζ)            # markov chain object 
    grid_ζ                    = mc_ζ.state_values                 # states 
    Π_ζ                       = mc_ζ.p                            # transition matrix 

    # other
    r                         = 0.04                              # exogeneous interest rate
end

# Results structure
mutable struct Results
    value_function::Array{Float64}        # value function
    policy_function::Array{Float64}       # savings policy function
end

function Initialize()
    @unpack N, N_a, N_ζ, N_ε = Primitives()

    # value_function, policy_function, and μ are 4-dimensional arrays
    # 1st dim is age, 2nd dim is asset holding, 3rd dim is persistent income shock, 4rd dim is transitory shock
    value_function       = zeros(N, N_a, N_ζ, N_ε)
    policy_function      = zeros(N, N_a, N_ζ, N_ε)

    Results(value_function, policy_function)

end

# Utility function
function u(c::Float64, γ::Float64)
    if (c > 0)
        (c^(1-γ))/(1 - γ)
    else
        -Inf
    end
end

function Solve_HH_problem!(results::Results; progress::Bool = false)
    @unpack β, N, γ, r, κ = Primitives()
    @unpack N_ε, N_ζ, Π_ε, Π_ζ, grid_ε, grid_ζ = Primitives()
    @unpack N_a, grid_a, grid_srl_a = Primitives()

    # In last period, agents consume everything.
    for (i_ζ,ζ) = enumerate(grid_ζ), (i_ε, ε) = enumerate(grid_ζ)
        results.value_function[N, :, i_ζ, i_ε] = u.(grid_a .+ exp(κ[N] + ζ + ε), γ)
    end

    # Backward induction starting at period before last period
    for t = (N-1):-1:1
        if (progress)
            println(t)
        end

        # iterates over assets, persistent shocks, and transitory shocks today
        for (i_ζ, ζ) = enumerate(grid_ζ), (i_ε, ε) = enumerate(grid_ζ)

            choice_lower = 1 # exploits monotonicity of policy function

            for (i_a, a) = enumerate(grid_a)

                # budget for savings and consumption
                income = exp(κ[t] + ζ + ε)
                budget = (1 + r) * a + income

                # stops grid search when value starts to decrease
                val_previous = -Inf

                # iterates over assets tomorrow
                for (i_a_p, a_p) = enumerate(grid_a)

                    # calculates utility
                    val = u(budget - a_p, γ)
                    val += β * Π_ζ[i_ζ, :]' * results.value_function[t+1, i_a_p, :, :] * Π_ε[i_ε, :]

                    # if utility starts to decrease, we save values and move on
                    if val < val_previous
                        results.value_function[t, i_a, i_ζ, i_ε] = val_previous
                        results.policy_function[t, i_a, i_ζ, i_ε] = grid_a[i_a_p-1]
                        choice_lower = i_a_p - 1
                        break

                    # if we're at the top of the grid, we save and move on
                    elseif i_a_p == N_a
                        results.value_function[t, i_a, i_ζ, i_ε] = val
                        results.policy_function[t, i_a, i_ζ, i_ε] = a_p
                    end

                    # update val_previous to check if utility is starting to decrease.
                    val_previous = val
                end
            end
        end
    end
end