####################################################################################################################
# ECON 899A Computational Economics
# Problem Set 6

# Krusell and Smith (1998)
# Alex von Hafften
# October 20, 2021

# This file defines the parameters of the model
####################################################################################################################

using Parameters

# Model primitives
@with_kw struct Primitives

    # Preference parameters
    β::Float64                = 0.8                                    # Discount rate
    A::Float64                = 1/200.0                                # Coefficient on labor in preferences

    # Production parameters
    θ::Float64                = 0.64                                   # Output elasticity of labor in production
    s_grid::Array{Float64, 1} = [3.98e-4, 3.58, 6.82, 12.18, 18.79]    # Shock grid
    n_s::Int64                = length(s_grid)                         # Number of states
    F::Array{Float64, 2}      = [0.6598 0.2600 0.0416 0.0331 0.0055;  
                                 0.1997 0.7201 0.0420 0.0326 0.0056; 
                                 0.2000 0.2000 0.5555 0.0344 0.0101;   # Transition function
                                 0.2000 0.2000 0.2502 0.3397 0.0101; 
                                 0.2000 0.2000 0.2500 0.3400 0.0100] 
    ν::Array{Float64, 1}      = [0.37, 0.4631, 0.1102, 0.0504, 0.0063] # Invariant distribution
    c_f::Float64              = 10.0                                   # Fixed cost
    c_e::Float64              = 5.0                                    # Entry cost

    # Other
    tolerence_EC::Float64     = 1e-4                                  # Tolerence for convergence of the entry condition
    tolerence_LMC::Float64    = 1e-4                                  # Tolerence for convergence of the labor market clearing condition
end

# Results structure
mutable struct Results
    p::Float64             # Final good price
    N_d::Array{Float64, 1} # labor demand
    π::Array{Float64, 1}   # Firm profits
    x::Array{Int64, 1}     # Exit decision
    W::Array{Float64, 1}   # Firm franchise value
    μ::Array{Float64, 1}   # Firm distribution
    M::Float64             # New entrant mass
end

# Initialize results structure
function Initialize()
    @unpack n_s = Primitives()
    p   = 1
    N_d = [1.3e-9, 10, 60, 300, 1000]
    π   = zeros(n_s)
    x   = fill(1, n_s)
    W   = zeros(n_s)
    μ   = ones(n_s)/n_s
    M   = 1

    return Results(p, N_d, π, x, W, μ, M)
end

# Solve for optimal labor demand by a given firm
function labor_demand(p::Float64, s::Float64, θ::Float64)
    return max((p*s*θ)^(1/(1-θ)), 0)
end

# returns static profits
function static_profit(p::Float64, s::Float64, n::Float64, θ::Float64, c_f::Float64)
    p*s*n^θ - n - p * c_f
end

# Bellman for exit decisions
function Bellman(P::Primitives, R::Results)

    # Initialize next policy and value function
    next_x = fill(0, P.n_s)
    next_W = zeros(P.n_s)

    # iterate over states today
    for i_s = 1:P.n_s

        # exiting today gives you current profit and zeros for all future periods
        W_exit = R.π[i_s]

        # not exiting gives you current profit plus discounted future franchise value
        W_stay = R.π[i_s]
        for i_s_p = 1:P.n_s
            W_stay += P.β * P.F[i_s, i_s_p] * R.W[i_s_p]
        end

        if (W_stay >= W_exit)
            next_x[i_s] = 0
            next_W[i_s] = W_stay
        else
            next_x[i_s] = 1
            next_W[i_s] = W_exit
        end
    end

    return next_x, next_W
end

function Solve_firm_problem(R::Results)
    P = Primitives()

    # Solve static labor demand decision
    R.N_d = labor_demand.(R.p, P.s_grid, P.θ)
    R.π   = static_profit.(R.p, P.s_grid, R.N_d, P.θ, P.c_f)

    # Solve dynamic firm exit decision
    i = 1
    err = 100

    while err > 0
        println(i)
        next_x, next_W = Bellman(P, R)
        err = sum(abs.(R.x .- next_x)) + sum(abs.(R.W .- next_W))
        println(err)
        R.x   = next_x
        R.W   = next_W
        i += 1
    end
    R
end