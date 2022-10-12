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
    N_lz_c::Int64                  # number of grid points
    MC_lz_c::MarkovChain           # tauchen as MarkovChain
    grid_lz_c::Vector{Float64}     # grid
    Π_lz_c::Matrix{Float64}        # transition probabilities

    # fine productivity
    N_lz_f::Int64                  # number of grid points
    MC_lz_f::MarkovChain           # tauchen as MarkovChain
    grid_lz_f::Vector{Float64}     # grid
    Π_lz_f::Matrix{Float64}        # transition probabilities

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

    # net worth
    min_w::Float64              # minimum
    max_w::Float64              # maximum
    N_w::Int64                  # number of grid points
    grid_w::Vector{Float64}     # grid
end

# creates Grids based on In_Primitives
function Initialize_Grids(P::Primitives)

    # coarse productivity
    N_lz_c        = 15
    MC_lz_c       = tauchen(N_lz_c, P.ρ, P.σ_ε, 0, 4)
    grid_lz_c     = collect(MC_lz_c.state_values)
    π_lz_c        = MC_lz_c.p
    min_lz        = grid_lz_c[1]
    max_lz        = grid_lz_c[N_lz_c]
    # min_lz and max_lz should be about the following...
    # min_lz        = -4*P.σ_ε/ sqrt(1- P.ρ^2)
    # max_lz        =  4*P.σ_ε/ sqrt(1- P.ρ^2)

    # fine productivity
    N_lz_f        = 60
    MC_lz_f       = tauchen(N_lz_f, P.ρ, P.σ_ε, 0, 4)
    grid_lz_f     = collect(MC_lz_f.state_values)
    π_lz_f        = MC_lz_f.p
    # max_lz and min_lz should be the same as above

    # capital
    k_bar      = ((P.δ / exp(max_lz))/P.α)^(1/(P.α-1))
    grid_k     = k_bar.*(1 - P.δ).^(15:-0.5:0)
    N_k        = length(grid_k)
    min_k      = grid_k[1]
    max_k      = grid_k[N_k]

    # debt
    min_b  = -(1-P.τ_c_p) * k_bar ^ P.α / P.r
    max_b  = (1-P.τ_c_p) * k_bar ^ P.α / P.r
    N_b    = Int64(floor(N_k/2))
    grid_b = collect(range(min_b, max_b; length = N_b))

    # net worth grid - HW doesn't say much about this...
    # min_w  = compute_nw(min_k, max_b, exp(min_lz), exp(min_lz), P.r, P) # approximate lower end w min k and z and max amount of debt with interest rate r
    # max_w  = compute_nw(max_k, min_b, exp(max_lz), exp(max_lz), P.r, P) # approximate upper end w max k and z and min amount of debt with interest rate r
    temp_grid_w = zeros(N_k, N_b, N_lz_c, N_lz_c)
    for (i_k_p, k_p) = enumerate(grid_k), (i_b_p, b_p) = enumerate(grid_b), (i_lz, lz) = enumerate(grid_lz_c), (i_lz_p, lz_p) = enumerate(grid_lz_c)
        y_p = exp(lz_p) * k_p ^ P.α - P.δ * k_p - P.r * b_p
        w_p = y_p - T_C(y_p, P) + k_p - b_p
        temp_grid_w[i_k_p, i_b_p, i_lz, i_lz_p] = w_p
    end
    min_w  = minimum(temp_grid_w)
    max_w  = maximum(temp_grid_w)
    N_w    = 100 # just a guess...
    grid_w = collect(range(min_w, max_w; length = N_w))

    return Grids(min_lz, max_lz, 
                 N_lz_c, MC_lz_c, grid_lz_c, π_lz_c, 
                 N_lz_f, MC_lz_f, grid_lz_f, π_lz_f, 
                 min_k, max_k, N_k, grid_k, 
                 min_b, max_b, N_b, grid_b,
                 min_w, max_w, N_w, grid_w)
end

# structure for Results
mutable struct Results
    r_tilde::Array{Float64} # bond prices
    vf::Array{Float64}      # value function
    pf_b::Array{Float64}    # bond policy function 
    pf_k::Array{Float64}    # capital policy function 
    pf_d::Array{Int64}    # default policy function
    w_bar::Array{Float64}   # default threshold net worth 
end

function Initialize_Results(P::Primitives, G::Grids)
    # start bond rates at r
    # bond price dimensions are capital, bonds, productivity
    r_tilde = zeros(G.N_k, G.N_b, G.N_lz_c)
    r_tilde .+= P.r

    # vf dimensions are realized net worth and current productivity
    vf   = zeros(G.N_w, G.N_lz_c)
    pf_b = zeros(G.N_w, G.N_lz_c)
    pf_k = zeros(G.N_w, G.N_lz_c)
    pf_d = fill(0, G.N_w, G.N_lz_c)
    w_bar = zeros(G.N_lz_c)

    Results(r_tilde, vf, pf_b, pf_k, pf_d, w_bar)
end
