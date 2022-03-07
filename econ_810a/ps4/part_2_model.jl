# Alex von Hafften
# ECON 810: Advanced Macro
# PS 4 - Part 2 - Code for solving model
# Professor Carter Braxton

using Parameters, Distributions, StatsBase

cd("/Users/vonhafften/Documents/UW Madison/problem_sets/econ_810a/ps4/")

# Primitive structure
@with_kw struct Primitives
    β::Float64           = 0.99                # hh discount rate
    r_f::Float64         = 0.04                # risk-free rate
    σ::Float64           = 2.0                 # coefficient of relative risk aversion
    α::Float64           = 0.70                # human capital technology parameter
    T::Int64             = 30                  # lifespan
    R::Vector{Float64}   = (1.0019).^(0:(T-1)) # rental rates of labor
end

# grids structure
@with_kw struct Grids

    # human capital grid
    min_h::Float64            = 0.01                              # lower bound
    max_h::Float64            = 10.0                              # upper bound
    N_h::Int64                = 100                               # number of points
    grid_srl_h                = range(min_h, max_h; length = N_h) # step range
    grid_h::Array{Float64, 1} = collect(grid_srl_h)               # grid array

    # savings grid
    min_k::Float64            = 0.01                              # lower bound
    max_k::Float64            = 10.0                              # upper bound
    N_k::Int64                = 100                               # number of points
    grid_srl_k                = range(min_k, max_k; length = N_k) # step range
    grid_k::Array{Float64, 1} = collect(grid_srl_k)               # grid array

    # time grid
    min_s::Float64            = 0.0                               # lower bound
    max_s::Float64            = 1.0                               # upper bound
    N_s::Int64                = 11                                # number of points
    grid_srl_s                = range(min_s, max_s; length = N_s) # step range
    grid_s::Array{Float64, 1} = collect(grid_srl_s)               # grid array

end

# Results
mutable struct Results
    # log normal human capital shocks distribution
    μ_z::Float64               # mean 
    σ_z::Float64               # sd 
    distribution_z             # distribution

    # initial human capital distribution 
    μ_h_0::Float64             # mean 
    σ_h_0::Float64             # sd 
    distribution_h_0           # distribution         

    vf::Array{Float64}         # value function
    pf_s::Array{Float64}       # value function
    pf_k::Array{Float64}       # value function
    
end

# Initialize results structure
function Initialize(μ_z, σ_z, μ_h_0, σ_h_0)
    @unpack T = Primitives()
    @unpack N_h, N_k = Grids()

    distribution_z = Normal(μ_z, σ_z)
    distribution_h_0 = Normal(μ_h_0, σ_h_0)

    # 1st dimension is age, 2nd dimension is human capital, 3rd dimension is assets.
    vf   = zeros(T, N_h, N_k)
    pf_s = zeros(T, N_h, N_k)
    pf_k = zeros(T, N_h, N_k)

    Results(μ_z, σ_z, distribution_z, μ_h_0, σ_h_0, distribution_h_0, vf, pf_s, pf_k)
end

# utility function
function u(c::Float64, σ::Float64)
    if (c > 0)
        return (c^(1-σ))/(1-σ)
    else
        return -1/eps()
    end
end

# human capital technology
function H(h::Float64, s::Float64, α::Float64)
    return h + (h*s)^α
end

# in terminal period all agents work and do not save
function Solve_terminal_period!(R::Results)
    P = Primitives()
    G = Grids()

    println(P.T)

    for (i_h, h) = enumerate(G.grid_h), (i_k, k) = enumerate(G.grid_k)
        
        consumption = P.R[P.T]*h + k*(1+P.r_f)
        R.vf[P.T, i_h, i_k] = u(consumption, P.σ)

    end
end

# solve nonterminal periods
function Solve_nonterminal_period!(R::Results)
    P = Primitives()
    G = Grids()

    # loop over periods
    for t = (P.T-1):-1:1
        println(t)

        # loop over state variables
        for (i_h, h) = enumerate(G.grid_h), (i_k, k) = enumerate(G.grid_k)
        
            candidate_max = -1/eps()

            # loop over choices
            for (i_s, s) = enumerate(G.grid_s), (i_k_p, k_p) = enumerate(G.grid_k)

                # consumption
                consumption = P.R[t] * h * (1-s) + k * (1 - P.r_f) - k_p

                if consumption <= 0
                    continue
                end

                # instantanteous value
                value = u(consumption, P.σ)

                # invert law of motion of human capital to get shocks
                z_p = log.(G.grid_h/H(h, s, P.α))
                z_p_breaks = midpoints(z_p)

                # probability of shocks
                z_p_Π = diff(vcat(0, cdf.(R.distribution_z, z_p_breaks), 1))

                # continuation value
                value += P.β * R.vf[t+1, :, i_k_p]' * z_p_Π

                if value >= candidate_max
                    R.vf[t, i_h, i_k] = value
                    R.pf_s[t, i_h, i_k] = s
                    R.pf_k[t, i_h, i_k] = k_p
                    candidate_max = value
                end
            end
        end
    end
end


function Solve!(R::Results)
    Solve_terminal_period!(R)
    Solve_nonterminal_period!(R)
end