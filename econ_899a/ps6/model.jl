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

using Parameters, LinearAlgebra, CSV, Tables, DataFrames

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
    c_e::Float64              = 5.0                                    # Entry cost

    # Other
    tolerence_EC::Float64     = 1e-4                                  # Tolerence for convergence of the entry condition
    tolerence_LMC::Float64    = 1e-4                                  # Tolerence for convergence of the labor market clearing condition
end

# Results structure
mutable struct Results
    # Parameters that change between simulations
    c_f::Float64           # Fixed cost
    α::Float64             # Type 1 Extreme Value distribution parameter 

    # First stage results
    p::Float64             # Final good price
    N_d::Array{Float64, 1} # labor demand
    π::Array{Float64, 1}   # Firm profits
    x::Array{Float64, 1}   # Exit decision
    W::Array{Float64, 1}   # Firm franchise value

    # Second stage results
    M::Float64             # Mass of entrants
    μ::Array{Float64, 1}   # Stationary firm distribution
    Π::Float64             # Aggregate profits
    L_d::Float64           # Aggregate labor demand
    L_s::Float64           # Aggregate labor supply
end

# Initialize results structure
function Initialize(; c_f::Float64 = 10.0, α::Float64 = -1.0)
    @unpack n_s = Primitives()

    price = 1.0
    N_d = zeros(n_s)
    π   = zeros(n_s)
    x   = fill(1, n_s)
    W   = zeros(n_s)
    M   = 1.0
    μ   = ones(n_s)
    Π   = 0.0
    L_d = 0.0
    L_s = 0.0

    return Results(c_f, α, price, N_d, π, x, W, M, μ, Π, L_d, L_s)
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
    next_x = zeros(P.n_s)
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
            next_x[i_s] = 0.0
            next_W[i_s] = W_stay

        else # exits if better than staying
            next_x[i_s] = 1.0
            next_W[i_s] = W_exit
        end
    end

    return next_x, next_W
end

# Bellman for exit decisions with type 1 EV shocks
function Exit_Bellman_random(P::Primitives, R::Results)

    # Initialize next policy and value function
    next_x = zeros(P.n_s)
    next_W = zeros(P.n_s)

    # iterate over states today
    for i_s = 1:P.n_s

        # exiting today gives you current profit and zeros for all future periods
        V_exit = R.π[i_s]

        # not exiting gives you current profit plus discounted future franchise value
        V_stay = R.π[i_s]
        for i_s_p = 1:P.n_s
            V_stay += P.β * P.F[i_s, i_s_p] * R.W[i_s_p]
        end

        # Using log-sum-exp trick
        c = max(R.α * V_stay, R.α * V_exit)
        next_W[i_s] = MathConstants.eulergamma / R.α + 1/R.α * (c + log(exp(R.α * V_stay - c) + exp(R.α * V_exit - c)))

        next_x[i_s] = exp(R.α * V_exit - c) / (exp.(R.α * V_stay - c) + exp(R.α * V_exit - c))
    end

    return next_x, next_W
end

# Solve firm problem
function Solve_firm_problem(R::Results)
    P = Primitives()

    # Solve static labor demand decision
    R.N_d = compute_labor_demand.(R.p, P.s_grid, P.θ)
    R.π   = compute_static_profit.(R.p, P.s_grid, R.N_d, P.θ, R.c_f)

    # Solve dynamic firm exit decision
    i = 1
    err = 100

    while err > 0

        # println(i) # for debugging
        # println(err) # for debugging
        
        # compute updated policy and value function guesses
        if (R.α == -1) # standard
            next_x, next_W = Exit_Bellman(P, R)
        else # with type 1 ev shocks
            next_x, next_W = Exit_Bellman_random(P, R)
        end 
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
function compute_entry_condition(W::Array{Float64, 1}, price::Float64)
    @unpack ν, c_e = Primitives()

    # Compute entry condition
    sum(W .* ν) / price - c_e
end

# Solves for price using bisection method.
function Solve_price(R::Results)
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
        # println(p_mid) # for debugging

        R.p = p_mid

        Solve_firm_problem(R)

        EC  = compute_entry_condition(R.W, p_mid) # evaluate entry condition at midpoint price

        if EC < 0 # EC is less than zero, move up lower bound.
            p_low = p_mid
        else # EC is larger than zero, move down upper bound.
            p_high = p_mid
        end

        p_mid = (p_high + p_low)/2 # compute new midpoint
        i += 1
    end

    R.p = p_mid

    R
end

####################################################################################################################
############################################### Solve for new entrant mass #########################################
####################################################################################################################

function compute_μ(R::Results)
    @unpack F, ν, n_s = Primitives()
    Z = reshape(repeat(1 .- R.x, n_s), n_s, n_s)' .* F'
    R.μ = R.M * inv(I - Z) * Z * ν
    R
end

# Computes μ using T_star operator to verify matrix algebra
function T_star(R::Results, P::Primitives)
    μ_p = zeros(P.n_s)

    for i_s = 1:P.n_s
        for i_s_p = 1:P.n_s
            μ_p[i_s_p] += (1 - R.x[i_s]) * P.F[i_s, i_s_p] * R.μ[i_s]
            μ_p[i_s_p] += (1 - R.x[i_s]) * P.F[i_s, i_s_p] * R.M * P.ν[i_s]
        end
    end

    return μ_p
end

function compute_μ_T_star(R::Results)
    P = Primitives()

    err, i = 100, 1

    while err > P.tolerence_LMC # borrow LMC tolerence_LMC
        μ_p = T_star(R, P) 
        err = maximum(abs.(R.μ .- μ_p))
        # println("Iteration #", i) # for debugging
        # println("Error is ", err) # for debugging
        R.μ = μ_p
        i += 1
    end
    R
end

# compute labor demand based on firm decision
function compute_labor_demand(R::Results)
    @unpack ν = Primitives()
    R.L_d = sum(R.N_d .* R.μ) + R.M * sum(R.N_d .* ν)
end

# compute aggregate firm profits
function compute_aggregate_profit(R::Results)
    @unpack ν = Primitives()
    R.Π = sum(R.π .* R.μ) + R.M * sum(R.π .* ν)
end

# compute labor supply based on FOC of HH problem
function compute_labor_supply(R::Results)
    @unpack A = Primitives()
    R.L_s = 1/A - compute_aggregate_profit(R)
end

# compute LMC
function compute_LMC(R::Results)
    compute_μ(R)
    return compute_labor_demand(R) - compute_labor_supply(R)
end 

# Solves for mass of new entrants using bisection method.
function Solve_M(R::Results)
    @unpack tolerence_LMC = Primitives()

    # Initial bounds and midpoint for m
    M_low = tolerence_LMC
    M_high = 10.0
    R.M = (M_high + M_low)/2

    # Loop variables
    LMC = 100
    i = 1

    # iterate until convergence
    while abs(LMC) > tolerence_LMC
        
        # println(i) # for debugging
        # println(LMC) # for debugging

        LMC  = compute_LMC(R) # evaluate entry condition at midpoint m

        if LMC < 0 # LMC is less than zero, move up lower bound.
            M_low = R.M
        else # LMC is larger than zero, move down upper bound.
            M_high = R.M
        end

        R.M = (M_high + M_low)/2 # compute new midpoint
        i += 1
    end
    R.M
end

# wrapper
function Solve_model(; c_f::Float64 = 10.0, α::Float64 = -1.0)
    R = Initialize(;c_f = c_f, α = α)
    Solve_price(R)
    Solve_M(R)
    R
end

################################################################################
############################ Functions to create summary table #################
################################################################################

process_results = function(results::Results)

    @unpack ν = Primitives()

    # rows
    c_f = results.c_f
    alpha = results.α
    price_level = results.p
    mass_incumbents = sum((1 .- results.x) .* results.μ)
    mass_entrants = results.M
    mass_exits = sum(results.x .* results.μ)
    aggregate_labor = results.L_d
    labor_incumbents = sum(results.N_d .* results.μ)
    labor_entrants = results.M * sum(results.N_d .* ν)
    frac_labor_entrants = results.M * sum(results.N_d .* ν)/ results.L_d

    # create vector of summary statistics
    [c_f, alpha, price_level, mass_incumbents, mass_entrants, mass_exits, 
    aggregate_labor, labor_incumbents, labor_entrants, frac_labor_entrants]
end

function create_table(results_vector::Array{Results})
    table = DataFrames.DataFrame(Tables.table(reduce(hcat,process_results.(results_vector))'))
    rename!(table, [:c_f, :alpha, :price_level, :mass_incumbents, :mass_entrants, :mass_exits, 
                    :aggregate_labor, :labor_incumbents, :labor_entrants, :frac_labor_entrants])
end