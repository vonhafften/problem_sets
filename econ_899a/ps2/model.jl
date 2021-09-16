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
using Parameters, Interpolations, Optim

# Primitive structure
@with_kw struct Primitives
    β::Float64                = 0.9932                           # Discount rate
    α::Float64                = 1.5                              # Coefficient of relative risk aversion
    S::Array{Float64, 1}      = [1.0, 0.5]                       # Earning states
    S_length::Int64           = length(S)                        # Number of earning states
    Π::Array{Float64, 2}      = [[0.97, 0.5] [0.03, 0.5]]        # Employment transition matrix
    a_min::Int64              = -2.0                             # Lower bound of savings grid
    a_max::Int64              = 5.0                              # Upper bound of savings grid
    a_length::Int64           = 1000                             # Number of points on asset grid
    a_grid_srl                = range(a_min, a_max;length = a_length)      # Savings grid step range
    a_grid::Array{Float64, 1} = collect(a_grid_srl)               # Savings grid array
end

mutable struct Results
    value_function::Array{Float64, 2}  # Value function`
    policy_function::Array{Float64, 2} # Policy function
    μ::Array{Float64, 1}               # Asset distribution
    q::Float64                         # Bond price
end

function Initialize()
    @unpack a_length, β = Primitives()

    value_function  = [zeros(a_length) zeros(a_length)]
    policy_function = [zeros(a_length) zeros(a_length)]
    μ               = ones(a_length)/a_length           # Start with uniform wealth distribution
    q               = (β+1)/2                           # We've assumed 1 > q > β, so start at midpoint

    Results(value_function, policy_function, μ, q)
end

################################################################################
######################### Functions to solve HH problem ########################
################################################################################

function HH_Bellman(results::Results)
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

            println(s, ", ", a)

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

function Solve_HH_problem(results::Results, tolerence::Float64 = 1e-4)

    err, i = 100.0, 0

    while err > tolerence
        print(i)
        i += 1
        v_next = HH_Bellman(results)
        err = abs.(maximum(v_next.-results.value_function))/abs(maximum(v_next))
        results.value_function = v_next
    end

    print("HH Problem converged in ", i, " iterations")
end

################################################################################
######################### Functions to find invariant asset distribution #######
################################################################################



################################################################################
######################### Functions to apply market clearing ###################
################################################################################



################################################################################
######################### Functions to run whole model #########################
################################################################################

function Solve_model(tolerence::Float64 = 1e-3)

    err, i = 100.0, 1

    while err > tolerence
        results = Initialize()

        Solve_HH_problem()
        # Solve HH bellman
        # Solve invariant asset distribution
        # Apply market clearing
    end

end
