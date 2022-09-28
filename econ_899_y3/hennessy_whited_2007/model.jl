# Code to solve Hennessy and Whited 2007

# Alex von Hafften
# Sept. 27

using Parameters, QuantEcon

# structure for model primitives that are externally calibrated
@with_kw struct Ex_Parameters
    # tax parameters
    τ_d_bar::Float64  = 0.12 # tax rate on cash distributions
    τ_i::Float64      = 0.29 # tax rate on interest income
    τ_c_p::Float64    = 0.40 # corporate tax rate on positive income
    τ_c_n::Float64    = 0.20 # corporate tax rate on negative income

    # other
    r::Float64        = 0.025 # risk-free rate 
    δ::Float64        = 0.15  # depreciation
end

# structure for model primitives that are internally calibrated
mutable struct In_Parameters
    α::Float64   # curvature of production function
    λ_0::Float64 # fixed cost of external equity issuance
    λ_1::Float64 # linear cost coefficient of -//-
    λ_2::Float64 # quadratic cost coefficient of -//-
    ξ::Float64   # bankruptcy cost parameter
    ϕ::Float64   # shape of distributions tax schedule
    σ_ε::Float64 # z shock variance
    ρ::Float64   # z drift 
end

# creates In_Parameters with estimated parameters from HW 2007
function Initialize_In_Parameters(;type="baseline")
    if type == "baseline"
        return In_Parameters(0.627, 0.598, 0.091, 0.0004, 0.104, 0.732, 0.118, 0.684)
    elseif type == "small"
        return In_Parameters(0.693, 0.951, 0.120, 0.0004, 0.151, 0.831, 0.159, 0.498)
    elseif type == "large"
        return In_Parameters(0.577, 0.389, 0.053, 0.0002, 0.084, 0.695, 0.086, 0.791)
    else
        error("Specify valid type in Initialize_In_Parameters()")
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
    π_lz_c::Matrix{Float64}        # transition probabilities

    # fine productivity
    N_lz_f::Int64                  # number of grid points
    MC_lz_f::MarkovChain           # tauchen as MarkovChain
    grid_lz_f::Vector{Float64}     # grid
    π_lz_f::Matrix{Float64}        # transition probabilities

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

# operating profits
function op(k::Float64, α::Float64)
    k^α
end

# inverse of operating profits
function op_inv(profits::Float64, α::Float64)
    profits^(1/α)
end

# operating profits derivative
function op_d(k::Float64, α::Float64)
    α * k ^ (α - 1)
end

# inverse of operating profits derivative
function op_d_inv(profits::Float64, α::Float64)
    (profits/α)^(1/(α-1))
end

# creates Grids based on In_Primitives
function Initialize_Grids(IP::In_Parameters)
    EP = Ex_Parameters()

    # coarse productivity
    N_lz_c        = 15
    MC_lz_c       = tauchen(N_lz_c, IP.ρ, IP.σ_ε, 0, 4)
    grid_lz_c     = collect(MC_lz_c.state_values)
    π_lz_c        = MC_lz_c.p
    min_lz        = grid_lz_c[1]
    max_lz        = grid_lz_c[N_lz_c]
    # min_lz and max_lz should be about the following...
    # min_lz        = -4*IP.σ_ε/ sqrt(1- IP.ρ^2)
    # max_lz        =  4*IP.σ_ε/ sqrt(1- IP.ρ^2)

    # fine productivity
    N_lz_f        = 60
    MC_lz_f       = tauchen(N_lz_f, IP.ρ, IP.σ_ε, 0, 4)
    grid_lz_f     = collect(MC_lz_f.state_values)
    π_lz_f        = MC_lz_f.p
    # max_lz and min_lz should be the same as above

    # capital
    k_bar      = op_d_inv(EP.δ/exp(max_lz), IP.α)
    grid_k     = k_bar.*(1-EP.δ).^(15:-0.5:0)
    N_k        = length(grid_k)
    min_k      = grid_k[1]
    max_k      = grid_k[N_k]

    # debt
    min_b  = -(1-EP.τ_c_p) * k_bar ^ IP.α / EP.r
    max_b  = (1-EP.τ_c_p) * k_bar ^ IP.α / EP.r
    N_b    = Int64(floor(N_k/2))
    grid_b = collect(range(min_b, max_b; length = N_b))

    return Grids(min_lz, max_lz, 
                 N_lz_c, MC_lz_c, grid_lz_c, π_lz_c, 
                 N_lz_f, MC_lz_f, grid_lz_f, π_lz_f, 
                 min_k, max_k, N_k, grid_k, 
                 min_b, max_b, N_b, grid_b)
end

