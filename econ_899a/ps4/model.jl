# Conesa and Krueger (1999)
# Alex von Hafften
# October 6, 2021

# ECON 899A Computational Economics
# Problem Set 4

# This file contains the functions to estimate the model.
# ./transition.jl calls this file and computes steady states of the model.
# ./steady_state.jl calls this file and computes transition paths of the model.

################################################################################
######################### Housekeeping functions ###############################
################################################################################

# Load libraries
using Parameters

# Primitive structure
@with_kw struct Primitives
    # life-cycle
    N::Int64                  = 66                               # lifespan
    Jᴿ::Int64                 = 46                               # retirement age
    n::Float64                = 0.011                            # population growth rate

    # preference parameters
    σ::Float64                = 2.0                              # coefficient of relative risk aversion
    γ::Float64                = 0.42                              # utility weight on consumption
    β::Float64                = 0.97                             # Discount rate

    # productivity
    η::Array{Float64, 1}      = map(x->parse(Float64,x), readlines("ef.txt")) # deterministic age-efficiency profile
    z::Array{Float64, 1}      = [3.0, 0.5]                       # Idiosyncratic productivity levels (must be only two states)
    e::Array{Float64, 2}      = η * z'                           # productivity levels
    z_length::Int64           = 2                                # Number of productivity states
    Π::Array{Float64, 2}      = [0.9261 0.0739; 0.0189 0.9811]   # Productivity persistance probabilities
    Π₀::Array{Float64, 1}     = [0.2037, 0.7963]                 # Erodic distribution of Π

    # production
    α::Float64                = 0.36                             # capital elasticity in production
    δ::Float64                = 0.06                             # capital depreciation rate

    # asset grid
    a_min::Float64            = 0.0                              # Lower bound of savings grid
    a_max::Float64            = 75.0                             # Upper bound of savings grid
    a_length::Int64           = 5000                             # Number of points on asset grid
    a_grid_srl                = range(a_min, a_max; length = a_length) # Savings grid step range
    a_grid::Array{Float64, 1} = collect(a_grid_srl)               # Savings grid array
end

# Utility function of retired agent
function uᴿ(c::Float64, γ::Float64, σ::Float64)
    if (c > 0)
        (c^((1-σ) * γ))/(1 - σ)
    else
        -1/eps()
    end
end

# labor supply decision
function labor_decision(a::Float64, a_p::Float64, e::Float64, θ::Float64,
                        γ::Float64, r::Float64, w::Float64)
    interior_solution = (γ*(1-θ)*e*w - (1-γ)*((1+r)*a - a_p)) / ((1-θ)*w*e)
    min(1, max(0, interior_solution))
end

# utility function of workers
function uᵂ(c::Float64, l::Float64, γ::Float64, σ::Float64)
    if (c > 0 && l >= 0 && l <= 1)
        (((c^γ) * ((1 - l)^(1-γ)))^(1-σ))/(1 - σ)
    else
        -1/eps()
    end
end
