####################################################################################################################

# FIN 971 Corporate Finance
# Problem Set 5

# Gomes (2001)
# Alex von Hafften
# December 15, 2021

# This file contains the model code.

####################################################################################################################

####################################################################################################################
################################## Define primitives and results structures ########################################
####################################################################################################################

using Parameters, QuantEcon

# Model parameters
@with_kw struct Primitives
    α_k::Float64              = 0.3                                    # Share of capital
    α_l::Float64              = 0.65                                   # Share of labor
    δ::Float64                = 0.145                                  # Depreciation rate
    β::Float64                = 1/1.065                                # Discount factor
    H::Float64                = 0.6                                    # Preference to leisure
    A::Float64                = 1.0                                    # Aggregate productivity
    f_p::Float64              = 0.01                                   # Fixed cost of production
end

@with_kw struct Grids
    # Productivity grid
    ρ::Float64                  = 0.762                                  # Persistency of productivity
    σ::Float64                  = 0.0352                                 # Std dev of productivity Shock
    z_n::Int64                  = 10                                     # Number of productivity states
    z_mc                        = tauchen(10, ρ, σ, 0, 4)                # productivity markov chain
    z_Π::Array{Float64, 2}      = z_mc.p                                 # productivity transition matrix
    z_grid::Array{Float64, 1}   = z_mc.state_values                      # productivity grid
    ϕ::Array{Float64, 1}        = stationary_distributions(z_mc)[1]      # stationary distribution of productivity process

    # Capital grid
    k_min::Float64              = 0.0                                    # Lower bound of capital grid
    k_max::Float64              = 5.0                                   # Upper bound of capital grid
    k_n::Int64                  = 51                                     # Number of capital grid points
    k_grid                      = range(k_min, k_max; length = k_n)      # Capital grid

end

# Results structure
mutable struct Results
    # Mutable parameters
    λ_0::Float64            # Floatation cost of external finance
    λ_1::Float64            # Variational cost of external finance

    # First stage results
    wage::Float64           # Final good price
    N_d::Array{Float64, 2}  # Firm-level labor demand
    π::Array{Float64, 2}    # Firm static profits
    x_pf::Array{Float64, 2} # Exit decision policy function
    k_pf::Array{Float64, 2} # Capital policy function
    vf::Array{Float64, 2}   # Firm franchise value

    # Second stage results
    B::Float64             # Mass of entrants
    μ::Array{Float64, 2}   # Stationary firm distribution
    Π::Float64             # Aggregate profits
    L_d::Float64           # Aggregate labor demand
    L_s::Float64           # Aggregate labor supply

    ###########################################################################
    # Other results
    ###########################################################################

    # matrix results

    cdf::Array{Float64, 2}  # CDF of Stationary firm distribution
    d_pf::Array{Float64, 2} # Dividend policy function
    i_pf::Array{Float64, 2} # Investment policy function
    λ::Array{Float64, 2}    # External financing cost

    breakdown_constrained::Array{Float64, 1} # Investment policy function
    
    Y::Float64              # aggregate output
    I::Float64              # aggregate investment
    Λ::Float64              # aggregate financing costs
    production_cost::Float64# aggregate production costs
    floatation_cost::Float64  # aggregate floatation costs

    size_mean::Float64      # average size
    i_k_mean::Float64       # mean of i/k
    i_k_std::Float64        # std of i/k
    tobin_q::Float64        # Mean q
    cash_flow_mean::Float64 # Mean cash flow
    cash_flow_std::Float64  # Std cash flow
    neg_i_frac::Float64     # Fraction of firms with negative investment
end

# Initialize results structure
function Initialize(λ_0::Float64, λ_1::Float64)
    @unpack k_n, z_n = Grids()

    # Mutable parameters

    # First stage results
    wage = 1.0
    N_d  = zeros(k_n, z_n)
    π    = zeros(k_n, z_n)
    x_pf = fill(1, k_n, z_n)
    k_pf = zeros(k_n, z_n)
    vf   = zeros(k_n, z_n)

    # Second stage results
    B    = 1.0
    μ    = ones(k_n, z_n)./(k_n*z_n)
    Π    = 0.0
    L_d  = 0.0
    L_s  = 0.0

    # Other results
    cdf  = zeros(k_n, z_n)
    d_pf = zeros(k_n, z_n)
    i_pf = zeros(k_n, z_n)
    λ    = zeros(k_n, z_n)
    breakdown_constrained = [0.0, 0.0, 0.0]

    return Results(λ_0, λ_1, wage, N_d, π, x_pf, k_pf, vf, B, μ, Π, L_d, L_s, cdf, d_pf, i_pf, λ, breakdown_constrained, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
end

####################################################################################################################
############################################### Solve firm problem #################################################
####################################################################################################################

# Solve for optimal labor demand by a given firm
function compute_labor_demand(wage::Float64, k::Float64, z::Float64, P::Primitives)
    return max((wage/(P.A * exp(z) * P.α_l * k^P.α_k))^(1/(P.α_l - 1)), 0)
end

# returns static profits
function compute_static_profit(wage::Float64, l::Float64, k::Float64, z::Float64, P::Primitives)
    return P.A * exp(z) * k^P.α_k * l^P.α_l - wage * l - P.f_p
end

# Bellman for exit decisions
function Bellman(P::Primitives, G::Grids, R::Results)

    # Initialize next policy and value function
    next_x_pf = zeros(G.k_n, G.z_n)
    next_k_pf = zeros(G.k_n, G.z_n)
    next_vf = zeros(G.k_n, G.z_n)

    # iterate over states today
    for (i_k, k) = enumerate(G.k_grid), (i_z, z) = enumerate(G.z_grid)
        candidate_max = -1/eps()

        # iterate over states tomorrow
        for (i_k_p, k_p) = enumerate(G.k_grid)

            # compute static values
            cash_flow = R.π[i_k, i_z]
            investment = (k_p - (1 - P.δ) * k)
            if (cash_flow > investment)
                funding_cost = 0
            else 
                funding_cost = R.λ_0 + R.λ_1 * (investment - cash_flow)
            end

            # compute value of exiting vs staying
            v_exit = cash_flow - investment - funding_cost + P.β * (1-P.δ) * k_p
            v_stay = cash_flow - investment - funding_cost
            for (i_z_p, z_p) = enumerate(G.z_grid)
                v_stay += P.β * G.z_Π[i_z, i_z_p] * R.vf[i_k_p, i_z_p]
            end

            # update if higher than candidate max
            if (v_exit > candidate_max)
                next_x_pf[i_k, i_z] = 0
                next_k_pf[i_k, i_z] = k_p
                next_vf[i_k, i_z] = v_exit
                candidate_max = v_exit
            end
            if (v_stay > candidate_max)
                next_x_pf[i_k, i_z] = 1
                next_k_pf[i_k, i_z] = k_p
                next_vf[i_k, i_z] = v_stay
                candidate_max = v_stay
            end
        end
    end
    return next_x_pf, next_k_pf, next_vf
end

# Solve firm problem
function Solve_firm_problem(R::Results)
    P = Primitives()
    G = Grids()

    # Solve static labor demand decision
    for i_k = 1:G.k_n, i_z = 1:G.z_n
        R.N_d[i_k, i_z] = compute_labor_demand(R.wage, G.k_grid[i_k], G.z_grid[i_z], P)
        R.π[i_k,   i_z] = compute_static_profit(R.wage, R.N_d[i_k, i_z], G.k_grid[i_k], G.z_grid[i_z], P)
    end

    # Solve dynamic firm exit decision
    i, err, maxiter = 1, 100, 100

    while (err > 1e-4) & (i < maxiter)

        # println(i) # for debugging
        # println(err) # for debugging
        
        # compute updated policy and value function guesses
        next_x_pf, next_k_pf, next_vf = Bellman(P, G, R)

        # sup norm
        err = maximum(vcat(abs.(R.x_pf .- next_x_pf), abs.(R.k_pf .- next_k_pf), abs.(R.vf .- next_vf)))
        
        # update
        R.x_pf   = next_x_pf
        R.k_pf   = next_k_pf
        R.vf     = next_vf
        i += 1
    end
    return R
end

####################################################################################################################
############################################### Solve for price ####################################################
####################################################################################################################

# Solves for price using bisection method.
function Solve_wage(R::Results)
    G = Grids()
    
    # Initial bounds and midpoint for price search using bisection method.
    wage_low = 0.75 
    wage_high = 1.5
    R.wage = (wage_high + wage_low)/2

    # Loop variables
    entry_condition, i, maxiter = 100, 1, 100

    # iterate until convergence
    while (abs(entry_condition) > 1e-3) & (i < maxiter)
        # println("Iteration #", i) # for debugging
        # println("The entry condition equals ", entry_condition) # for debugging
        # println("Current wage guess is ", R.wage) # for debugging

        # evaluate entry condition at midpoint wage
        R = Solve_firm_problem(R)
        entry_condition  = G.ϕ' * R.vf[1,:]

        # EC is less than zero, move up lower bound.
        if entry_condition < 0 
            wage_high = R.wage
        else # EC is larger than zero, move down upper bound.
            wage_low = R.wage
        end

        R.wage = (wage_high + wage_low)/2 # compute new midpoint
        i += 1
    end

    return R
end

####################################################################################################################
############################################### Solve for new entrant mass #########################################
####################################################################################################################

# Solves stationary distribution for B = 1
function compute_μ(R::Results)
    P = Primitives()
    G = Grids()

    transition_matrix = zeros(G.k_n * G.z_n, G.k_n * G.z_n)

    for (i_k, k) in enumerate(G.k_grid), (i_z, z) in enumerate(G.z_grid)
        for (i_k_p, k_p) in enumerate(G.k_grid), (i_z_p, z_p) in enumerate(G.z_grid)

            if k_p == R.k_pf[i_k, i_z]
                row = i_k + G.k_n * (i_z - 1)
                col = i_k_p + G.k_n * (i_z_p - 1)

                transition_matrix[row, col] += R.x_pf[i_k, i_z] * G.z_Π[i_z, i_z_p] 
            end
        end
    end

    transition_matrix_entrant = zeros(G.z_n, G.k_n * G.z_n)

    for (i_z, z) in enumerate(G.z_grid)
        for (i_k_p, k_p) in enumerate(G.k_grid), (i_z_p, z_p) in enumerate(G.z_grid)
            if k_p == R.k_pf[1, i_z]
                row = i_z
                col = i_k_p + G.k_n * (i_z_p - 1)

                transition_matrix_entrant[row, col] += R.x_pf[1, i_z] * G.z_Π[i_z, i_z_p] 
            end
        end
    end    

    μ = reshape(R.μ, G.k_n * G.z_n, 1)
    μ_p = zeros(G.k_n * G.z_n, 1)

    err = 100.0

    while err > 0.001
        # again we assume B = 1
        μ_p = transition_matrix' * μ + transition_matrix_entrant' * G.ϕ

        err = maximum(abs.(μ_p .- μ))
        μ = μ_p
    end

    return reshape(μ, G.k_n, G.z_n)
end

# Solves for mass of new entrants using bisection method.
function Solve_B(R::Results)
    @unpack H = Primitives()
    @unpack ϕ = Grids()

    μ = compute_μ(R)

    labor_demand_1 = sum(R.π .* R.x_pf .* μ) + sum(R.π[1,:] .* ϕ)
    aggregate_profit_1 = sum(R.N_d .* μ .* R.x_pf) + sum(R.N_d[1,:] .* ϕ)
    
    R.B = 1/(H*(labor_demand_1 + aggregate_profit_1/R.wage))
    R.μ = μ * R.B
    R.L_d = sum(R.π .* R.x_pf .* R.μ) + R.B * sum(R.π[1,:] .* ϕ)
    R.Π = sum(R.N_d .* R.μ .* R.x_pf) + R.B * sum(R.N_d[1,:] .* ϕ)
    R.L_s = 1/H - R.Π/R.wage
    
    if abs(R.L_s - R.L_d) > 1e-5
        println("Labor supply: ", R.L_s)
        println("Labor demand: ", R.L_d)
        error("Labor market clearing not working.")
    end

    return R

end

# wrapper
function Solve_model(λ_0::Float64, λ_1::Float64)
    P = Primitives()
    G = Grids()

    R = Initialize(λ_0, λ_1)

    print("Solving for wage... ")
    Solve_wage(R)
    println("Done.")
    
    print("Solving for entrant mass... ")
    Solve_B(R)
    println("Done.")

    # Other results
    k_matrix = G.k_grid * ones(G.z_n)'
    z_matrix = ones(G.k_n) * G.z_grid' 

    # matrix results
    R.cdf = cumsum(R.μ; dims = 1)./sum(R.μ; dims = 1)
    R.i_pf = R.k_pf .- (1 - P.δ) .* k_matrix
    R.λ = (R.i_pf .- R.π .> 0) .* (λ_0 .+ λ_1 .* (R.i_pf .- R.π))
    R.d_pf = R.π .- R.i_pf .- R.λ

    # aggregate moments
    R.breakdown_constrained[1] = sum((R.d_pf .< -0.1) .* R.μ)/sum(R.μ)
    R.breakdown_constrained[2] = sum((R.d_pf .> -0.1) .* (R.d_pf .< 0.1) .* R.μ)/sum(R.μ)
    R.breakdown_constrained[3] = sum((R.d_pf .>  0.1) .* R.μ)/sum(R.μ)
    R.Y = sum((P.A .* exp.(z_matrix) .* k_matrix .^ P.α_k .* R.N_d .^ P.α_l .- P.f_p) .* R.μ .* R.x_pf) - R.B * P.f_p
    R.I = sum(R.i_pf .* R.x_pf .* R.μ) + R.B * sum(R.k_pf[1,:] .* G.ϕ)
    R.Λ = sum(R.λ .* R.x_pf .* R.μ) + R.B * sum(R.λ[1,:] .* G.ϕ)
    R.production_cost = sum(R.μ .* P.f_p)
    R.floatation_cost = sum((R.λ .> 0) .* λ_0 .* R.μ)

    # cross-sectional moments
    M = sum(R.μ)

    R.size_mean = 1/(M - R.B) * sum(k_matrix .* R.μ .* R.x_pf)
    R.i_k_mean = 1/(M - R.B) * sum((k_matrix .> 0) .* (R.i_pf ./ k_matrix) .* R.μ .* R.x_pf)
    R.i_k_std = sqrt(1/(M - R.B) * sum((k_matrix .> 0) .* ((R.i_pf ./ k_matrix) .- R.i_k_mean).^2 .* R.μ .* R.x_pf))
    R.tobin_q = 1/(M - R.B) * sum((k_matrix .> 0.0) .* (R.vf ./ k_matrix) .* R.μ .* R.x_pf)
    R.cash_flow_mean = 1/(M - R.B) * sum((k_matrix .> 0.0) .* (R.π ./ k_matrix) .* R.μ .* R.x_pf)
    R.cash_flow_std = sqrt(1/(M - R.B) * sum((k_matrix .> 0) .* ((R.π ./ k_matrix) .- R.cash_flow_mean).^2 .* R.μ .* R.x_pf))
    R.neg_i_frac = sum((R.i_pf .< 0) .* R.μ .* R.x_pf)

    return R
end




