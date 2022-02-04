# Kaplan and Violante (2010)
# Alex von Hafften
# January 31, 2022

# ECON 810A Advanced Macro Theory
# Problem set 1 - Part 2

# This file contains the functions to simulate the model.
# ./part_2_run.jl calls this file and runs the model.

# Load libraries
using Parameters, Statistics

include("part_2_model.jl")

# Results structure
mutable struct Simulation
    trials::Int64               # Number of individuals simulated
    ζ_states::Array{Int64}      # persistent income shock states
    ε_states::Array{Int64}      # transitory income shock states
    assets::Array{Float64}      # Assets
    income::Array{Float64}      # income
    consumption::Array{Float64} # consumption
end

# Initializes the simulation structure
function Initialize_Simulation(trials::Int64)
    @unpack N = Primitives()

    ζ_states = zeros(trials, N)
    ε_states = zeros(trials, N)
    assets = zeros(trials, N)
    income = zeros(trials, N)
    consumption = zeros(trials, N)

    Simulation(trials, ζ_states, ε_states, assets, income, consumption)

end

# simulates the markov chains
function simulate_shocks!(simulation_results::Simulation)
    @unpack N, mc_ε, mc_ζ = Primitives()

    # simulate markov chains
    for i = 1:simulation_results.trials
        simulation_results.ζ_states[i,:] = simulate_indices(mc_ζ, N)
        simulation_results.ε_states[i,:] = simulate_indices(mc_ε, N)
    end

end

function compute_assets!(simulation_results::Simulation, results::Results)
    @unpack N, grid_a = Primitives()
    @unpack grid_ζ, grid_ε, κ = Primitives()
    @unpack r = Primitives()

    for i = 1:simulation_results.trials
        for t = 1:(N-1)

            # find indices
            i_a = argmin(abs.(simulation_results.assets[i,t] .- grid_a))
            i_ζ = simulation_results.ζ_states[i, t]
            i_ε = simulation_results.ε_states[i, t]

            # apply policy function to find in tomorrows assets
            simulation_results.assets[i,t+1] = results.policy_function[t, i_a, i_ζ, i_ε]

            # fill in income
            simulation_results.income[i,t] = exp(κ[t] + grid_ζ[i_ζ] + grid_ε[i_ε])

            # fill in consumption
            assets_today = simulation_results.assets[i,t] * (1 + r)
            assets_tomorrow = simulation_results.assets[i,t+1]
            income = simulation_results.income[i,t]

            # I'm getting negative consumption...
            simulation_results.consumption[i,t] = max(assets_today + income - assets_tomorrow, eps())
        end

        # terminal period
        i_ζ = simulation_results.ζ_states[i, N]
        i_ε = simulation_results.ε_states[i, N]

        # get income shocks
        ζ = grid_ζ[i_ζ]
        ε = grid_ε[i_ε]

        # fill in income
        simulation_results.income[i,N] = exp(κ[N] + ζ + ε)

        # fill in consumption
        assets_today = simulation_results.assets[i,N] * (1 + r)
        income = simulation_results.income[i,N]
        simulation_results.consumption[i,N] = assets_today + income
    end
end


# computes the bpp shock variance
function compute_bbp_coefficients(simulation::Simulation)
    @unpack N, κ, ρ_ζ = Primitives()

    # create vectors for pseudo difference (plus lags and leads)
    delta_y = zeros((N-3)*simulation.trials)
    delta_y_lead = zeros((N-3)*simulation.trials)
    delta_y_lag = zeros((N-3)*simulation.trials)
    delta_c = zeros((N-3)*simulation.trials)

    for i = 1:simulation.trials

        # pull out life-cycle earning component
        y_i = log.(simulation.income[i,:]) - κ
        c_i = simulation.consumption[i,:]

        # compute difference in consumption and pseudo difference in income
        delta_y_i = y_i[2:end] - (ρ_ζ * y_i[1:(end-1)])
        delta_c_i = log.(c_i[2:end]) - log.(c_i[1:(end-1)])

        # indices to store this observations changes.
        indices = (((i-1)*(N-3))+1):((i*(N-3)))

        # store pseudo differences income
        delta_y_lag[indices]  = delta_y_i[1:(end-2)]
        delta_y[indices]      = delta_y_i[2:(end-1)]
        delta_y_lead[indices] = delta_y_i[3:end] 

        # store consumption differences
        delta_c[indices] = delta_c_i[2:(end-1)]
    end

    # Using formulas from slide 80 from Carter's presentation
    cov_y_y_lead = cov(delta_y, delta_y_lead)
    cov_c_y_lead = cov(delta_c, delta_y_lead)
    cov_y_y_lead_lag = cov(delta_y, ρ_ζ^2 .* delta_y_lag .+ ρ_ζ .* delta_y .+ delta_y_lead)
    cov_c_y_lead_lag = cov(delta_c, ρ_ζ^2 .* delta_y_lag .+ ρ_ζ .* delta_y .+ delta_y_lead)
    
    var_ε = (-1/ρ_ζ) * cov_y_y_lead
    var_ζ = (1/ρ_ζ) * cov_y_y_lead_lag

    pass_ε = cov_c_y_lead/cov_y_y_lead
    pass_ζ = cov_c_y_lead_lag/cov_y_y_lead_lag

    return (var_ε, var_ζ, pass_ε, pass_ζ)
end