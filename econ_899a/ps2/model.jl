# Huggett (1993)
# Alex von Hafften
# September 16, 2021

# ECON 899A Computational Economics
# Problem set 2

# This file contains the functions to estimate the model.
# ./run.jl calls this file and runs the model.

################################################################################
######################### Housekeeping functions ###############################
################################################################################

# Load libraries
using Parameters, Interpolations, Optim, LinearAlgebra

# Primitive structure
@with_kw struct Primitives
    β::Float64                = 0.9821                           # Discount rate
    α::Float64                = 2.5                              # Coefficient of relative risk aversion
    S::Array{Float64, 1}      = [0.03, 0.01]                       # Earning states
    S_length::Int64           = length(S)                        # Number of earning states
    Π::Array{Float64, 2}      = [[0.95, 0.05] [0.95, 0.05]]        # Employment transition matrix
    a_min::Float64              = -0.3                             # Lower bound of savings grid
    a_max::Float64              = 5.0                              # Upper bound of savings grid
    a_length::Int64           = 1001                             # Number of points on asset grid
    a_grid_srl                = range(a_min, a_max;length = a_length)      # Savings grid step range
    a_grid::Array{Float64, 1} = collect(a_grid_srl)               # Savings grid array
end

mutable struct Results
    value_function::Array{Float64, 2}  # Value function`
    policy_function::Array{Float64, 2} # Policy function
    μ::Array{Float64, 2}               # Bond holdings distribution
    q::Float64                         # Bond price
end

function Initialize()
    @unpack a_length, β = Primitives()

    value_function  = [zeros(a_length) zeros(a_length)]
    policy_function = [zeros(a_length) zeros(a_length)]
    μ               = [ones(a_length) ones(a_length)] / (a_length*2)    # Start with uniform bond holding distribution
    q               = (β+1)/2                                           # We've assumed 1 > q > β, so start at midpoint

    Results(value_function, policy_function, μ, q)
end

################################################################################
######################### Functions to solve HH problem ########################
################################################################################

function HH_Bellman(results::Results; progress::Bool = false)
    @unpack a_min, a_max, a_length, a_grid, a_grid_srl, S, S_length, Π, α, β = Primitives()

    v_next = zeros(a_length, S_length)  # Initialize next value function iteration

    # Interpolate value function
    v_interp_e = scale(interpolate(results.value_function[:,1], BSpline(Cubic(Line(OnGrid())))), a_grid_srl)
    v_interp_u = scale(interpolate(results.value_function[:,2], BSpline(Cubic(Line(OnGrid())))), a_grid_srl)

    # cycle over employment statuses
    for i_s = 1:S_length

        s = S[i_s]

        # loop over asset holdings
        for i_a = 1:a_length

            a = a_grid[i_a] # value of assets

            if progress
                println(s, ", ", a)
            end

            budget = s + a # budget

            val(a_p) = ((budget - a_p * results.q)^(1-α)-1)/(1-α) + β * Π[i_s, 1] * v_interp_e(a_p) + β * Π[i_s, 2] * v_interp_u(a_p)
            obj(a_p) = - val(a_p)

            opt = optimize(obj, a_min, min(budget, a_max)) # Optimize using Brent's method

            results.policy_function[i_a, i_s] = opt.minimizer
            v_next[i_a, i_s] = -opt.minimum
        end
    end
    v_next
end

function Solve_HH_problem(results::Results, tolerence::Float64 = 1e-5; progress::Bool = false)

    err, i = 100.0, 0

    while err > tolerence
        i += 1

        if progress
            println(i)
        end

        v_next = HH_Bellman(results)
        err = abs.(maximum(v_next.-results.value_function))/abs(maximum(v_next))
        results.value_function = v_next
    end

    println("HH Problem converged in ", i, " iterations")
end

################################################################################
################## Functions to find invariant bond holding distribution #######
################################################################################

function Apply_policy_function_to_μ(results::Results; progress::Bool = false)
    @unpack a_length, S_length, Π = Primitives()

    μ_next = zeros(a_length, S_length)

    for i_s = 1:S_length
        for i_a = 1:a_length

            if progress
                println(i_s, ", ", i_a)
            end

            a_p = results.policy_function[i_a, i_s]

            i_a_p_e = argmin(abs.(a_p .- a_grid))
            i_a_p_u = argmin(abs.(a_p .- a_grid))

            μ_next[i_a_p_e, 1] = μ_next[i_a_p_e, 1] + results.μ[i_a, i_s] * Π[i_s, 1]
            μ_next[i_a_p_u, 2] = μ_next[i_a_p_u, 2] + results.μ[i_a, i_s] * Π[i_s, 2]

        end
    end
    μ_next
end

function Solve_invariant_μ(results::Results, tolerence::Float64 = 1e-5; progress::Bool = false)

    err, i = 100.0, 0

    while err > tolerence

        i += 1

        if progress
            println(i)
        end

        μ_next = Apply_policy_function_to_μ(results)
        err = abs.(maximum(μ_next .- results.μ))/abs(maximum(μ_next))
        results.μ = μ_next
    end

    println("Invariant μ converged in ", i, " iterations")

end

################################################################################
################ Function to update price based on market clearing  ############
################################################################################

function Update_price(results::Results, tolerence::Float64 = 1e-3)
    @unpack a_grid, β = Primitives()

    excess_demand = -sum(results.μ .* [a_grid a_grid])

    if excess_demand > tolerence
        q_hat = results.q + (β - results.q)/2*abs(excess_demand)

        println("Excess Demand is positive: ", excess_demand)
        println("Lower bond price from ", results.q, " to ", q_hat)

        results.q = q_hat

        return(false)
    elseif excess_demand < -tolerence
        q_hat = results.q + (1 - results.q)/2*abs(excess_demand)

        println("Excess Demand is negative: ", excess_demand)
        println("Raise bond price from ", results.q, " to ", q_hat)

        results.q = q_hat

        return(false)
    else
        println("Excess Demand is within tolerence: ", excess_demand)

        return(true)
    end
end

################################################################################
######################### Functions to run whole model #########################
################################################################################

function Solve_model(results::Results)

    converged = false

    while !converged
        Solve_HH_problem(results)
        Solve_invariant_μ(results)
        converged = Update_price(results)
    end
end

################################################################################
######################### Function calculate wealth distribution ###############
################################################################################

function calculate_wealth_distribution(results::Results)
    @unpack a_grid, S, S_length, a_length = Primitives()

    w = zeros(a_length, S_length)

    for i_s = 1:S_length
        for i_a = 1:a_length

            i_w = argmin(abs.(a_grid[i_a] .+ S[i_s] .- a_grid))

            w[i_w, i_s] = results.μ[i_a, i_s]
        end
    end
    w
end

function calculate_lorenz_curve(w::Array{Float64, 2})
    @unpack a_grid, a_length = Primitives()

    x = cumsum(w[:,1] .+ w[:,2])
    y = cumsum((w[:,1] .+ w[:,2]) .* a_grid)

    unique([x/x[a_length] y/y[a_length]]; dims = 1)
end

# https://en.wikipedia.org/wiki/Gini_coefficient
function calculate_gini(l::Array{Float64, 2})
    widths = diff(l[:,1])
    heights = ((l[1:end-1,1] .+ l[2:end,1])./2 .- (l[1:end-1,2] .+ l[2:end,2])./2)
    a = sum(widths .* heights)

    l_pos = l[l[:,2].>0, :]
    widths = diff(l_pos[:,1])
    heights = (l_pos[1:end-1,2] .+ l_pos[2:end,2])./2
    b = sum(widths .* heights)

    a/(a+b)
end

################################################################################
######################### Function calculate welfare change ####################
################################################################################

function calculate_w_fb()
    @unpack α, β, S, Π = Primitives()
    stationary_distribution = (Π^1000000)[1, :]
    c_fb = stationary_distribution[1] * S[1] + stationary_distribution[2] * S[2]
    ((c_fb)^(1 - α) - 1)/((1 - α) * (1 - β))
end

function calculate_λ(results::Results, w_fb::Float64)
    @unpack α, β = Primitives()

    numerator = w_fb + 1 /((1 - α)*(1 - β))
    denominator = results.value_function .+ (1 ./((1 .- α).*(1 .- β)))

    (numerator./denominator).^(1/(1 .- α)) .- 1
end
