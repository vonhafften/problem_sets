####################################################################################################################
# ECON 899A Computational Economics
# Problem Set 6

# Hopenhayn and Rogerson (1993)
# Alex von Hafften
# October 20, 2021

# This file contains the model code
####################################################################################################################

####################################################################################################################
################################## Define primitives and results structures ########################################
####################################################################################################################

using Parameters, LinearAlgebra

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
end

# Initialize results structure
function Initialize(price::Float64)
    @unpack n_s = Primitives()
    N_d = [1.3e-9, 10, 60, 300, 1000]
    π   = zeros(n_s)
    x   = fill(1, n_s)
    W   = zeros(n_s)

    return Results(price, N_d, π, x, W)
end

####################################################################################################################
############################################### Solve firm problem #################################################
####################################################################################################################

# Solve for optimal labor demand by a given firm
function compute_labor_demand(p::Float64, s::Float64, θ::Float64)
    return max((p*s*θ)^(1/(1-θ)), 0)
end

# returns static profits
function compute_static_profit(p::Float64, s::Float64, n::Float64, θ::Float64, c_f::Float64)
    p*s*n^θ - n - p * c_f
end

# Bellman for exit decisions
function Exit_Bellman(P::Primitives, R::Results)

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

        # stays if better than exiting
        if (W_stay >= W_exit)
            next_x[i_s] = 0
            next_W[i_s] = W_stay

        else # exits if better than staying
            next_x[i_s] = 1
            next_W[i_s] = W_exit
        end
    end

    return next_x, next_W
end

# Solve firm problem
function Solve_firm_problem(R::Results)
    P = Primitives()

    # Solve static labor demand decision
    R.N_d = compute_labor_demand.(R.p, P.s_grid, P.θ)
    R.π   = compute_static_profit.(R.p, P.s_grid, R.N_d, P.θ, P.c_f)

    # Solve dynamic firm exit decision
    i = 1
    err = 100

    while err > 0

        # println(i) # for debugging
        # println(err) # for debugging
        
        # compute updated policy and value function guesses
        next_x, next_W = Exit_Bellman(P, R)
        
        # sup norm
        err = sum(abs.(R.x .- next_x)) + sum(abs.(R.W .- next_W))
        
        # update
        R.x   = next_x
        R.W   = next_W
        i += 1
    end
    R
end

####################################################################################################################
############################################### Solve for price ####################################################
####################################################################################################################

# Given a price, this function computes the entry condition
function compute_entry_condition(price::Float64)
    @unpack ν, c_e = Primitives()

    # Solve firm problem at price
    R = Solve_firm_problem(Initialize(price))

    # Compute entry condition
    sum(R.W .* ν) / R.p - c_e
end

# Solves for price using bisection method.
function Solve_price()
    @unpack tolerence_EC = Primitives()

    # Initial bounds and midpoint for price search using bisection method.
    p_low = tolerence_EC
    p_high = 10.0
    p_mid = (p_high + p_low)/2

    # Loop variables
    EC = 100
    i = 1

    # iterate until convergence
    while abs(EC) > tolerence_EC
        # println(i) # for debugging
        # println(EC) # for debugging

        EC  = compute_entry_condition(p_mid) # evaluate entry condition at midpoint price

        if EC < 0 # EC is less than zero, move up lower bound.
            p_low = p_mid
        else # EC is larger than zero, move down upper bound.
            p_high = p_mid
        end

        p_mid = (p_high + p_low)/2 # compute new midpoint
        i += 1
    end

    p_mid
end

####################################################################################################################
############################################### Solve for new entrant mass #########################################
####################################################################################################################

# Computes μ using T_star operator to verify other function
function T_star(μ::Array{Float64, 1}, R::Results, M::Float64, P::Primitives)
    μ_p = zeros(P.n_s)

    for i_s = 1:P.n_s
        for i_s_p = 1:P.n_s
            μ_p[i_s_p] += (1 - R.x[i_s]) * P.F[i_s, i_s_p] * μ[i_s]
            μ_p[i_s_p] += (1 - R.x[i_s]) * P.F[i_s, i_s_p] * M * P.ν[i_s]
        end
    end

    return μ_p
end

function compute_μ(R::Results, M::Float64)
    P = Primitives()

    err, i = 100, 1

    μ = ones(P.n_s)

    while err > P.tolerence_LMC # borrow LMC tolerence_LMC
        μ_p = T_star(μ, R, M, P) 
        err = maximum(abs.(μ .- μ_p))
        # println("Iteration #", i) # for debugging
        # println("Error is ", err) # for debugging
        μ = μ_p
        i += 1
    end
    μ
end

# compute labor demand based on firm decision
function compute_labor_demand(R::Results, M::Float64, μ::Array{Float64, 1})
    @unpack ν = Primitives()
    sum(R.N_d .* μ) + M * sum(R.N_d .* ν)
end

# compute labor supply based on FOC of HH problem
function compute_labor_supply(R::Results, M::Float64, μ::Array{Float64, 1})
    @unpack A, ν = Primitives()
    Π = sum(R.π .* μ) + M * sum(R.π .* ν)
    return 1/A - Π
end

# compute LMC
function compute_LMC(R::Results, M::Float64)
    μ = compute_μ(R, M)
    return compute_labor_demand(R, M, μ) - compute_labor_supply(R, M, μ)
end 

# Solves for mass of new entrants using bisection method.
function Solve_M(R::Results)
    @unpack tolerence_LMC = Primitives()

    # Initial bounds and midpoint for m
    m_low = tolerence_LMC
    m_high = 10.0
    m_mid = (m_high + m_low)/2

    # Loop variables
    LMC = 100
    i = 1

    # iterate until convergence
    while abs(LMC) > tolerence_LMC
        
        # println(i) # for debugging
        # println(LMC) # for debugging

        LMC  = compute_LMC(R, m_mid) # evaluate entry condition at midpoint m

        if LMC < 0 # LMC is less than zero, move up lower bound.
            m_low = m_mid
        else # LMC is larger than zero, move down upper bound.
            m_high = m_mid
        end

        m_mid = (m_high + m_low)/2 # compute new midpoint
        i += 1
    end

    m_mid
end