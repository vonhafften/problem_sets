# Code to solve Hennessy and Whited 2007

# Alex von Hafften
# Sept. 27, 2022

using Parameters, QuantEcon

# structure for model primitives that are externally calibrated
@with_kw struct Ex_Primitives

    # tax parameters
    τ_d_bar::Float64  = 0.12 # tax rate on cash distributions
    τ_i::Float64      = 0.29 # tax rate on individual income
    τ_c_p::Float64    = 0.40 # tax rate on positive corporate income
    τ_c_n::Float64    = 0.20 # tax rate on negative corporate income

    # other
    r::Float64        = 0.025 # risk-free rate 
    δ::Float64        = 0.15  # depreciation
end

# structure for model primitives that are externally and internally calibrated
mutable struct Primitives

    # tax parameters
    τ_d_bar::Float64 # tax rate on cash distributions
    τ_i::Float64     # tax rate on interest income
    τ_c_p::Float64   # tax rate on positive corporate income
    τ_c_n::Float64   # tax rate on negative corporate income

    # other
    r::Float64 # risk-free rate 
    δ::Float64 # depreciation

    # internally calibrated primitives
    α::Float64   # curvature of production function
    λ_0::Float64 # fixed cost of external equity issuance
    λ_1::Float64 # linear cost coefficient of -//-
    λ_2::Float64 # quadratic cost coefficient of -//-
    ξ::Float64   # bankruptcy cost parameter
    ϕ::Float64   # shape of distributions tax schedule
    σ_ε::Float64 # z shock variance
    ρ::Float64   # z drift 
end

# creates Primitives with estimated parameters from HW 2007
function Initialize_Primitives(;type="baseline")
    EP = Ex_Primitives()
    if type == "baseline"
        return Primitives(EP.τ_d_bar, EP.τ_i, EP.τ_c_p, EP.τ_c_n, EP.r, EP.δ, 0.627, 0.598, 0.091, 0.0004, 0.104, 0.732, 0.118, 0.684)
    elseif type == "small"
        return Primitives(EP.τ_d_bar, EP.τ_i, EP.τ_c_p, EP.τ_c_n, EP.r, EP.δ, 0.693, 0.951, 0.120, 0.0004, 0.151, 0.831, 0.159, 0.498)
    elseif type == "large"
        return Primitives(EP.τ_d_bar, EP.τ_i, EP.τ_c_p, EP.τ_c_n, EP.r, EP.δ, 0.577, 0.389, 0.053, 0.0002, 0.084, 0.695, 0.086, 0.791)
    else
        error("Specify valid type in Initialize_Primitives()")
    end
end

# structure for Grids
mutable struct Grids
    # productivity
    min_lz::Float64                # minimum
    max_lz::Float64                # maximum

    # coarse productivity
    N_lz::Int64                  # number of grid points
    MC_lz::MarkovChain           # tauchen as MarkovChain
    grid_lz::Vector{Float64}     # grid
    Π_lz::Matrix{Float64}        # transition probabilities

    # capital
    min_k::Float64              # minimum
    max_k::Float64              # maximum
    N_k::Int64                  # number of grid points
    grid_k::Vector{Float64}     # grid

    # debt
    min_b::Float64              # minimum
    max_b::Float64              # maximum
    N_b::Int64                  # number of grid points
    grid_b::Vector{Float64}     # grid

end

# creates Grids based on In_Primitives
function Initialize_Grids(P::Primitives)

    # coarse productivity
    N_lz        = 15
    MC_lz       = tauchen(N_lz, P.ρ, P.σ_ε, 0, 4)
    grid_lz     = collect(MC_lz.state_values)
    Π_lz        = MC_lz.p
    min_lz      = grid_lz[1]    # should be about min_lz = -4*P.σ_ε/ sqrt(1- P.ρ^2)
    max_lz      = grid_lz[N_lz] # should be about max_lz =  4*P.σ_ε/ sqrt(1- P.ρ^2)    

    # capital
    k_bar      = ((P.δ / exp(max_lz))/P.α)^(1/(P.α-1))
    grid_k     = k_bar.*(1 - P.δ).^(15:-0.5:0)
    N_k        = length(grid_k)
    min_k      = grid_k[1]
    max_k      = grid_k[N_k]

    # debt
    N_b    = Int64(floor(N_k/2))
    max_b  = (1-P.τ_c_p) * k_bar ^ P.α / P.r
    min_b  = -(1-P.τ_c_p) * k_bar ^ P.α / P.r
    temp_grid = collect(range(0.001, 0.999; length = N_b))
    temp_grid = log.(temp_grid ./ (1 .- temp_grid))
    grid_b = temp_grid ./ temp_grid[1] * min_b

    return Grids(min_lz, max_lz, 
                 N_lz, MC_lz, grid_lz, Π_lz, 
                 min_k, max_k, N_k, grid_k, 
                 min_b, max_b, N_b, grid_b)
end

# structure for Results
mutable struct Results

    # net worth grid
    min_w::Float64              # minimum
    max_w::Float64              # maximum
    N_w::Int64                  # number of grid points
    grid_w::Vector{Float64}     # grid

    # results
    q::Array{Float64}       # bond prices
    vf::Array{Float64}      # value function
    pf_b::Array{Float64}    # bond policy function 
    pf_k::Array{Float64}    # capital policy function 
    w_bar::Array{Float64}   # default threshold net worth 
    lz_d::Array{Float64}      # default policy function
end

# net worth grid - HW doesn't say much about this...
function compute_w_grid(q::Array{Float64}, N_w::Int64, P::Primitives, G::Grids)
    temp_grid_w = zeros(G.N_k, G.N_b, G.N_lz, G.N_lz)
    
    for (i_k_p, k_p) = enumerate(G.grid_k), (i_b_p, b_p) = enumerate(G.grid_b), (i_lz, lz) = enumerate(G.grid_lz), (i_lz_p, lz_p) = enumerate(G.grid_lz)
        y_p = exp(lz_p) * k_p ^ P.α - P.δ * k_p - (1-q[i_k_p,i_b_p,i_lz]) * b_p
        w_p = y_p - T_C(y_p, P) + k_p - q[i_k_p,i_b_p,i_lz]*b_p
        temp_grid_w[i_k_p, i_b_p, i_lz, i_lz_p] = w_p
    end

    min_w  = minimum(temp_grid_w) - 0.1
    max_w  = maximum(temp_grid_w) + 0.1
    
    linear_grid = collect(range(0.0, log(max_w-min_w + 1); length = N_w))

    return exp.(linear_grid) .+ min_w .- 1
end

function Initialize_Results(P::Primitives, G::Grids)
    
    # start bond rates at r
    # bond price dimensions are capital, bonds, productivity
    q = zeros(G.N_k, G.N_b, G.N_lz)
    q .+= 1/(1+P.r)

    #  w grid
    N_w    = 20 # just a guess...
    grid_w = compute_w_grid(q, N_w, P, G)
    min_w  = minimum(grid_w)
    max_w  = maximum(grid_w)

    # vf dimensions are realized net worth and current productivity
    vf   = zeros(N_w, G.N_lz)
    pf_b = zeros(N_w, G.N_lz)
    pf_k = zeros(N_w, G.N_lz)
    w_bar = zeros(G.N_lz)
    lz_d = zeros(G.N_k, G.N_b, G.N_lz)

    Results(min_w, max_w, N_w, grid_w, q, vf, pf_b, pf_k, w_bar, lz_d)
end
