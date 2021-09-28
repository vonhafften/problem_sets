# Conesa and Krueger (1999)
# Alex von Hafften
# September 22, 2021

# ECON 899A Computational Economics
# Problem set 3

# This file contains the functions to estimate the model.
# ./run.jl calls this file and runs the model.

################################################################################
######################### Housekeeping functions ###############################
################################################################################

# Load libraries
using Parameters, Interpolations, Optim

# Primitive structure
@with_kw struct Primitives
    N::Int64                  = 66                               # lifespan
    n::Float64                = 0.011                            # population growth rate
    a₁::Float64               = 0.0                              # initial assets
    Jᴿ::Int64                 = 46                               # retirement age
    θ::Float64                = 0.11                             # proportional labor tax
    γ::Float64                = 0.42                             # weight on consumption in utility function
    σ::Float64                = 2.0                              # coefficient of relative risk aversion
    η::Array{Float64, 1}      = map(x->parse(Float64,x), readlines("ef.txt")) # deterministic age-efficiency profile
    z::Array{Float64, 1}      = [3.0, 0.5]                       # Idiosyncratic productivity levels
    z_length::Int64           = length(z)                        # Number of productivity states
    Π::Array{Float64, 2}      = [0.9261 0.0739; 0.0189 0.9811]   # Productivity persistance probabilities
    Π₀::Array{Float64, 1}     = [0.2037, 0.7963]                 # Erodic distribution of Π
    e::Array{Float64, 2}      = η * z'                           # productivity levels
    α::Float64                = 0.36                             # capital elasticity in production
    δ::Float64                = 0.06                             # capital depreciation rate
    β::Float64                = 0.97                             # Discount rate
    a_min::Float64            = 0.01                             # Lower bound of savings grid
    a_max::Float64            = 75.0                              # Upper bound of savings grid
    a_length::Int64           = 5000                             # Number of points on asset grid
    a_grid_srl                = range(a_min, a_max; length = a_length) # Savings grid step range
    a_grid::Array{Float64, 1} = collect(a_grid_srl)               # Savings grid array
end

mutable struct Results
    value_function::Array{Float64}        # value function
    policy_function::Array{Float64}       # savings policy function
    labor_supply::Array{Float64}          # Optimal labor supply
    w::Float64                            # wage
    r::Float64                            # interest rate
    b::Float64                            # pension benefit
end

function Initialize()
    @unpack a_length, z_length, N, Jᴿ = Primitives()

    value_function       = zeros(N, a_length, z_length)
    policy_function      = zeros(N, a_length, z_length)
    labor_supply         = zeros(Jᴿ, a_length, z_length)
    w = 1.05
    r = 0.05
    b = 0.2

    Results(value_function, policy_function, labor_supply, w, r, b)
end

################################################################################
######################### Functions to solve HH problem ########################
################################################################################

function uᴿ(c::Float64, γ::Float64, σ::Float64)
    if (c > 0)
        (c^((1-σ) * γ))/(1 - σ)
    else
        -Inf
    end
end

function Solve_retiree_problem(results::Results; progress::Bool = false)
    @unpack β, N, Jᴿ, z_length, e, γ, σ = Primitives()
    @unpack a_min, a_max, a_length, a_grid, a_grid_srl = Primitives()

    results.value_function[N,:, 1] = uᴿ.(a_grid, γ, σ)

    for j = (N-1):-1:Jᴿ
        if (progress)
            println(j)
        end

        value_function_i = scale(interpolate(results.value_function[j+1,:, 1], BSpline(Cubic(Line(OnGrid())))), a_grid_srl)

        for i_a = 1:a_length

            budget = (1 + results.r) * a_grid[i_a] + results.b

            val(a_p) = uᴿ(budget - a_p, γ, σ) + β * value_function_i(a_p)
            obj(a_p) = - val(a_p)

            opt = optimize(obj, a_min, min(budget, a_max)) # Optimize using Brent's method

            results.policy_function[j, i_a, 1] = opt.minimizer
            results.value_function[j, i_a, 1] = -opt.minimum

        end
    end

    for i_z = 2:z_length
        results.policy_function[:, :, i_z] = results.policy_function[:, :, 1]
        results.value_function[:, :, i_z] = results.value_function[:, :, 1]
    end
end


function labor_decision(a::Float64, a_p::Float64, e::Float64, θ::Float64,  γ::Float64, r::Float64, w::Float64)
    interior_solution =  (γ * (1-θ) * e * w - (1- γ)*((1+r)*a - a_p)) / ((1- θ) * w * e)
    min(1, max(0, interior_solution))
end

function uᵂ(c::Float64, l::Float64, γ::Float64, σ::Float64)
    if (c > 0 && l >= 0 && l <= 1)
        (((c^γ) * ((1 - l)^(1-γ)))^(1-σ))/(1 - σ)
    else
        -Inf
    end
end

function Solve_worker_problem(results::Results; progress::Bool = false)
    @unpack β, Π, Jᴿ, θ, e, z_length, γ, σ = Primitives()
    @unpack a_min, a_max, a_length, a_grid, a_grid_srl = Primitives()

    for j = (Jᴿ-1):-1:1
        if (progress)
            println(j)
        end

        for i_z = 1:z_length
            choice_lower = 1 # exploits monotonicity of policy function
            for i_a = 1:a_length

                val_previous = -Inf # exploits monotonicity of value function

                for i_a_p = choice_lower:a_length

                    l = labor_decision(a_grid[i_a], a_grid[i_a_p], e[j, i_z], θ, γ, results.r, results.w)

                    budget = results.w * (1 - θ) * e[j, i_z] * l + (1 + results.r) * a_grid[i_a]
                    c = budget - a_grid[i_a_p]

                    val = uᵂ(c, l, γ, σ)

                    for i_z_p = 1:z_length
                        val += β * Π[i_z, i_z_p] * results.value_function[j+1, i_a_p, i_z_p]
                    end

                    if val < val_previous || i_a == a_length
                        results.value_function[j, i_a, i_z] = val_previous
                        results.policy_function[j, i_a, i_z] = a_grid[i_a_p-1]
                        results.labor_supply[j, i_a, i_z] = labor_decision(a_grid[i_a], a_grid[i_a_p-1], e[j, i_z], θ, γ, results.r, results.w)
                        choice_lower = i_a_p - 1
                        break
                    end
                    val_previous = val
                end
            end
        end
    end
end

function Solve_HH_problem(results::Results; progress::Bool = false)

    Solve_retiree_problem(results)
    Solve_worker_problem(results)

end
