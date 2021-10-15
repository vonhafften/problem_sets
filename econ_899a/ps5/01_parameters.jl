#############################################################################
# ECON 899A Computational Economics
# Problem set 5

# Krusell and Smith (1998)
# Alex von Hafften
# October 13, 2021

# Based on code written by Phil Coyle

# This file defines the parameters of the model including model primitives, grids, and shock structures.

#############################################################################

using Parameters

# Compute stationary distribution of any Markov process
function compute_Π_star(Π::Array{Float64})    
    x = length(Π[1,:])
    
    Π_star = zeros(x)
    Π_tmp  = Π^10000

    for i=1:x
        Π_star[i] = Π_tmp[i,i]
    end

    Π_star
end

# Model primitives
@with_kw struct Primitives
    # Preference parameters
    β::Float64                = 0.99                    # discount rate

    # Production technology parameters
    α::Float64                = 0.36                    # elasticity of capital in production
    δ::Float64                = 0.025                   # Capital depreciation
    e_bar::Float64            = 0.3271                  # Labor efficiency per unit of time worked
end


@with_kw struct Simulation
    # Simulation parameters
    T::Int64                  = 11000                   # Number of periods in simulation
    N::Int64                  = 5000                    # Number of individuals in simulation
    burn_in::Int64            = 1001                    # First period after burn-in period
    λ::Float64                = 0.5                     # Adjustment parameter
    tol_vfi::Float64          = 1e-4                    # Tolerence for the value function iteration
    tol_coef::Float64         = 1e-4                    # Tolerence for coefficient (a_0, etc.) iteration
    tol_r2::Float64           = 1.0 - 1e-2              # Tolerence for R^2 iteration
    max_iterations::Int64     = 10000                   # Maximum iterations
end

# Grids for aggregate productivity, employment, aggregate capital, and individual capital
@with_kw struct Grids

    # Aggregate productivity
    z::Array{Float64, 1}      = [1.01, 0.99]            # TFP states
    n_z::Int64                = length(z)               # number of states

    # Employment
    ε::Array{Float64, 1}      = [0, 1]                  # Idiosyncratic employment opportunities
    n_ε::Int64                = length(ε)               # number of states

    # Individual capital grid
    k_min::Float64                = 0.0                               # Minimum of individual capital grid
    k_max::Float64                = 15.0                              # Maximum of individual capital grid
    n_k::Int64                    = 61                                # Number of individual capital grid points
    k_grid_srl::Array{Float64, 1} = range(k_min, k_max; length = n_k) # individual capital grid step range length
    k_grid::Array{Float64, 1}     = collect(k_grid_srl)               # individual capital grid

    # Aggregate capital grid
    K_min::Float64                = 11.0                              # Minimum of aggregate capital grid
    K_max::Float64                = 15.0                              # Maximum of aggregate capital grid
    n_K::Int64                    = 17                                # Number of aggregate capital grid points
    K_grid_srl::Array{Float64, 1} = range(K_min, K_max; length = n_K) # Aggregate capital grid step range length
    K_grid::Array{Float64, 1}     = collect(K_grid_srl)               # Aggregate capital grid
end


# Transition matrices for business cycles and employment.
@with_kw struct Shocks
    # parameters of aggregate state transition matrix:
    d_g::Float64  = 8.0     # Duration (Good Times)
    d_b::Float64  = 8.0     # Duration (Bad Times)

    # transition probabilities for aggregate states
    pgg::Float64 = (d_g-1.0)/d_g
    pgb::Float64 = 1.0 - (d_b-1.0)/d_b
    pbg::Float64 = 1.0 - (d_g-1.0)/d_g
    pbb::Float64 = (d_b-1.0)/d_b

    # Transition matrix for aggregate state
    Π_z::Array{Float64,2} = [pgg pgb; pbg pbb]

    # Stationary aggregate state distribution
    Π_z_star::Array{Float64, 1} = compute_Π_star(Π_z)

    # parameters of aggregate state transition matrix:
    d_ug::Float64 = 1.5     # Unemp Duration (Good Times)
    d_ub::Float64 = 2.5     # Unemp Duration (Bad Times)
    u_g::Float64  = 0.04    # Fraction Unemp (Good Times)
    u_b::Float64  = 0.1     # Fraction Unemp (Bad Times)

    # transition probabilities for aggregate states and staying unemployed
    pgg00::Float64 = (d_ug-1.0)/d_ug
    pbb00::Float64 = (d_ub-1.0)/d_ub
    pbg00::Float64 = 1.25*pbb00
    pgb00::Float64 = 0.75*pgg00

    # transition probabilities for aggregate states and becoming employed
    pgg01::Float64 = (u_g - u_g*pgg00)/(1.0-u_g)
    pbb01::Float64 = (u_b - u_b*pbb00)/(1.0-u_b)
    pbg01::Float64 = (u_b - u_g*pbg00)/(1.0-u_g)
    pgb01::Float64 = (u_g - u_b*pgb00)/(1.0-u_b)

    # transition probabilities for aggregate states and becoming unemployed
    pgg10::Float64 = 1.0 - (d_ug-1.0)/d_ug
    pbb10::Float64 = 1.0 - (d_ub-1.0)/d_ub
    pbg10::Float64 = 1.0 - 1.25*pbb00
    pgb10::Float64 = 1.0 - 0.75*pgg00

    # transition probabilities for aggregate states and staying employed
    pgg11::Float64 = 1.0 - (u_g - u_g*pgg00)/(1.0-u_g)
    pbb11::Float64 = 1.0 - (u_b - u_b*pbb00)/(1.0-u_b)
    pbg11::Float64 = 1.0 - (u_b - u_g*pbg00)/(1.0-u_g)
    pgb11::Float64 = 1.0 - (u_g - u_b*pgb00)/(1.0-u_b)

    # Employment Transition Matrices conditional on aggregate state.
    Π_gg::Array{Float64,2} = [pgg11 pgg01; pgg10 pgg00]
    Π_bg::Array{Float64,2} = [pgb11 pgb01; pgb10 pgb00]
    Π_gb::Array{Float64,2} = [pbg11 pbg01; pbg10 pbg00]
    Π_bb::Array{Float64,2} = [pbb11 pbb01; pbb10 pbb00]

    # Uncoditional employment transition matrix
    Π_ε::Array{Float64,2} = [pgg*Π_gg pgb*Π_gb; pbg*Π_bg pbb*Π_bb]

    # Stationary employment distribution
    Π_ε_star::Array{Float64, 1} = compute_Π_star(Π_ε)
end

