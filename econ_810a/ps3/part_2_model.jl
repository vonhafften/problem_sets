# Alex von Hafften
# ECON 810: Advanced Macro
# PS 3 - Part 2 - Code for solving model
# Professor Carter Braxton

using Parameters

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_810a/ps3/")

# Primitive structure
@with_kw struct Primitives
    β_f::Float64 = 0.99  # firm discount rate
    β::Float64   = 0.99  # hh discount rate
    r_f::Float64 = 0.04  # risk-free rate
    δ::Float64   = 0.1   # job destruction rate
    ζ::Float64   = 1.6   # matching elasticity
    κ::Float64   = 0.995 # vacancy posting cost
    z::Float64   = 0.4   # unemployment benefit
    σ::Float64   = 2.0   # coefficient of relative risk aversion
    p_L::Float64 = 0.5   # probability of human capital drop while unemployed
    p_H::Float64 = 0.05  # probability of human capital improvement while employed
    T::Int64     = 120   # lifespan
    τ::Float64   = 0.2   # marginal tax rate
end

# grids structure
@with_kw struct Grids

    # human capital grid
    min_h::Float64            = 0.5                               # lower bound
    max_h::Float64            = 1.5                               # upper bound
    N_h::Int64                = 21                                # number of points
    grid_srl_h                = range(min_h, max_h; length = N_h) # step range
    grid_h::Array{Float64, 1} = collect(grid_srl_h)               # grid array

    # wage piece rate grid
    min_ω::Float64            = 0.0                               # lower bound
    max_ω::Float64            = 1.0                               # upper bound
    N_ω::Int64                = 21                               # number of points
    grid_srl_ω                = range(min_ω, max_ω; length = N_ω) # step range
    grid_ω::Array{Float64, 1} = collect(grid_srl_ω)               # grid array

    # bond grid
    min_b::Float64            = 0.0                               # lower bound
    max_b::Float64            = 30.0                              # upper bound
    N_b::Int64                = 51                               # number of points
    grid_srl_b                = range(min_b, max_b; length = N_b) # step range
    grid_b::Array{Float64, 1} = collect(grid_srl_b)               # grid array

end

# Results
mutable struct Results
    j_vf::Array{Float64}       # firm value function
    θ::Array{Float64}          # market tightness

    w_vf::Array{Float64}       # employed worker value function
    w_b_pf::Array{Float64}     # employed worker asset policy function

    u_vf::Array{Float64}       # unemployed worker value function
    u_b_pf::Array{Float64}     # unemployed worker asset policy function
    u_ω_pf::Array{Float64}     # unemployed worker piece rate policy function
end

# Initialize results structure
function Initialize()
    @unpack T = Primitives()
    @unpack N_ω, N_b, N_h = Grids()

    # 1st dimension is age, 2nd dimension is wage piece rate, 3rd dimension is human capital.
    j_vf = zeros(T, N_ω, N_h)
    θ    = zeros(T, N_ω, N_h)
    
    w_vf   = zeros(T, N_ω, N_b, N_h)
    w_b_pf = zeros(T, N_ω, N_b, N_h)

    u_vf   = zeros(T, N_b, N_h)
    u_b_pf = zeros(T, N_b, N_h)
    u_ω_pf = zeros(T, N_b, N_h)
    
    Results(j_vf, θ, w_vf, w_b_pf, u_vf, u_b_pf, u_ω_pf)
end

# utility function
function u(c::Float64, σ::Float64)
    if (c > 0)
        return (c^(1-σ))/(1-σ)
    else
        return -1/eps()
    end
end

# production function
function f(h::Float64)
    if (h > 0)
        return h
    else
        return 0
    end
end

# job finding rate
# input: market tightness 
# output: job finding rate 
function p(θ::Float64, ζ::Float64)
    return θ/((1 + θ^ζ)^(1/ζ))
end

# hiring rate
# input: market tightness
# output hiring rate
function p_f(θ::Float64, ζ::Float64)
    return 1/(1 + θ^(ζ))^(1/ζ)
end

# hiring rate inverse
# input: hiring rate
# output: market tightness
function p_f_inv(p_f::Float64, ζ::Float64)
    return (p_f^(-ζ) - 1)^(1/ζ)
end

# free entry to compute market tightness θ
# input: firm value and fixed cost of posting vacancy
# output: market tightness
function free_entry(j::Float64, κ::Float64, ζ::Float64)
    # nonbinding free entry condition 
    # (i.e. if even at hiring rate equal to one, but value is less than fixed cost of posting vacancy)
    if κ > j
        return 0.0
    
    # binding free entry condition
    else 
        return p_f_inv(κ/j, ζ)
    end
end

function Solve_terminal_period!(R::Results; progress::Bool = false)
    P = Primitives()
    G = Grids()

    # cycle over levels of human capital
    for (i_h, h) = enumerate(G.grid_h)
        # cycle over piece rates
        for (i_ω, ω) = enumerate(G.grid_ω)

            # compute firm value
            R.j_vf[P.T, i_ω, i_h] = (1-ω) * f(h)
            
            # invert free entry condition to compute market tightness
            R.θ[P.T, i_ω, i_h] = free_entry(R.j_vf[P.T, i_ω, i_h], P.κ, P.ζ)
            
            # compute employed worker value (eats everything)
            # cycle over asset holdings
            for (i_b, b) = enumerate(G.grid_b)
                R.w_vf[P.T, i_ω, i_b, i_h] = u((1-P.τ)*ω*f(h) + b, P.σ)
            end
        end

        # compute unemployed work value (eats everything)
        for (i_b, b) = enumerate(G.grid_b)
            R.u_vf[P.T, i_b, i_h] = u(P.z + b, P.σ)
        end
    end

    if progress
        println(P.T)
    end
end

function Solve_nonterminal_periods!(R::Results; progress::Bool = false)
    P = Primitives()
    G = Grids()

    # nonterminal periods
    for t = (P.T-1):-1:1

        # cycle over levels of human capital
        for (i_h, h) = enumerate(G.grid_h)
            
            # find index for human capital changes
            i_h_up = min(i_h + 1, G.N_h)
            i_h_down = max(i_h - 1, 1)
            
            # cycle over piece rates
            for (i_ω, ω) = enumerate(G.grid_ω)

                # compute firm value
                # instanteous value
                R.j_vf[t, i_ω, i_h] = (1-ω) * f(h)
                # continuation value if human capital increases
                R.j_vf[t, i_ω, i_h] += P.β_f * (1-P.δ) * P.p_H * R.j_vf[t+1, i_ω, i_h_up]
                # continuation value if human capital stays the same
                R.j_vf[t, i_ω, i_h] += P.β_f * (1-P.δ) * (1-P.p_H) * R.j_vf[t+1, i_ω, i_h]
            
                # invert free entry condition to compute market tightness
                R.θ[t, i_ω, i_h] = free_entry(R.j_vf[t, i_ω, i_h], P.κ, P.ζ)
            
                # compute employed worker value
                # cycle over asset holdings
                for (i_b, b) = enumerate(G.grid_b)
                    # compute budget
                    budget = (1-P.τ)*ω*f(h) + b
                    candidate_max = -1/eps()

                    # cycle over possible asset holdings tomorrow
                    for (i_b_p, b_p) = enumerate(G.grid_b)
                        # consumption is budget minus cost of bonds
                        consumption = budget - (1/(1+P.r_f)) *b_p

                        # instanteous value
                        value = u(consumption, P.σ)

                        # continuation value if stays employed and human capital stays same
                        value += P.β * (1-P.δ) * (1-P.p_H) * R.w_vf[t+1, i_ω, i_b_p, i_h]

                        # continuation value if stays employed and human capital increases
                        value += P.β * (1-P.δ) * P.p_H * R.w_vf[t+1, i_ω, i_b_p, i_h_up]

                        # continuation value if becomes unemployed and human capital stays same
                        value += P.β * P.δ * (1-P.p_H) * R.u_vf[t+1, i_b_p, i_h]
                        
                        # continuation value if becomes unemployed and human capital increases
                        value += P.β * P.δ * P.p_H * R.u_vf[t+1, i_b_p, i_h_up]

                        # save if better than candidate
                        if value > candidate_max
                            R.w_vf[t, i_ω, i_b, i_h] = value
                            R.w_b_pf[t, i_ω, i_b, i_h] = b_p
                            candidate_max = value
                        end
                    end
                end
            end
        
            # compute unemployed work value (eats everything)
            # cycle over assets today
            for (i_b, b) = enumerate(G.grid_b)
                # budget is savings plus unemployment benefit
                budget = P.z + b
                candidate_max = -1/eps()

                # cycle over assets tomorrow
                # cycle over piece rates
                for (i_b_p, b_p) = enumerate(G.grid_b), (i_ω_p, ω_p) = enumerate(G.grid_ω)
                    consumption = budget - (1/(1+P.r_f)) * b_p

                    # instanteous value
                    value = u(consumption, P.σ)

                    # job finding rates dependent on 
                    jfr_flat = p(R.θ[t+1, i_ω_p, i_h], P.ζ)
                    jfr_down = p(R.θ[t+1, i_ω_p, i_h_down], P.ζ)

                    # continuation value stays unemployed and human capital stays the same
                    value += P.β * (1-jfr_flat) * (1-P.p_L) * R.u_vf[t+1, i_b_p, i_h]
                    
                    # continuation value stays unemployed and human capital decreases
                    value += P.β * (1-jfr_down) * P.p_L * R.u_vf[t+1, i_b_p, i_h_down]

                    # continuation value becomes employed and human capital stays the same
                    value += P.β * jfr_flat * (1-P.p_L) * R.w_vf[t+1, i_ω_p, i_b_p, i_h]

                    # continuation value becomes employed and human capital decreases
                    value += P.β * jfr_down * P.p_L * R.w_vf[t+1, i_ω_p, i_b_p, i_h_down]

                    # save if better than candidate
                    if value > candidate_max
                        R.u_vf[t, i_b, i_h]  = value
                        R.u_b_pf[t, i_b, i_h] = b_p
                        R.u_ω_pf[t, i_b, i_h] = ω_p
                        candidate_max = value
                    end
                end
            end
        end
        if progress
            println(t)
        end
    end
end

function Solve!(R::Results; progress::Bool = false)
    Solve_terminal_period!(R; progress = progress)
    Solve_nonterminal_periods!(R; progress = progress)
end
