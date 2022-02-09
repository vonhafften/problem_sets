# Ljungqvist and Sargent (1998)
# Alex von Hafften
# January 8, 2022

# ECON 810A Advanced Macro Theory
# Problem Set 2 - Part 2

using Parameters, Distributions

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_810a/ps2/")

# Primitive structure
@with_kw struct Primitives

    # parameters
    T::Int64                  = 360                               # lifespan
    r::Float64                = 0.04                              # exogeneous interest rate
    β::Float64                = (1/(1+r))^(1/12)                  # Discount rate
    δ::Float64                = 0.033                             # Probability of being laid off
    ψ_e::Float64              = 0.5                               # probability of hc increase during employment
    ψ_u::Float64              = 0.2                               # probability of hc decrease during unemployment
    b::Float64                = 0.5                               # unemployment benefit

    # human capital grid
    min_h::Float64            = 1.0                               # Lower bound
    max_h::Float64            = 2.0                               # Upper bound
    N_h::Int64                = 25                                # Number of points
    Δ::Float64                = 1/25                              # human capital grid step size
    grid_srl_h                = range(min_h, max_h; length = N_h) # Step range
    grid_h::Array{Float64, 1} = collect(grid_srl_h)               # Grid array

    # search effort grid grid
    min_s::Float64            = 0.0                               # Lower bound
    max_s::Float64            = 1.0                               # Upper bound
    N_s::Int64                = 41                                # Number of points
    grid_srl_s                = range(min_s, max_s; length = N_s) # Step range
    grid_s::Array{Float64, 1} = collect(grid_srl_s)               # Grid array

    # wage offer grid
    μ_w::Float64              = 0.5                               # mean
    σ_w::Float64              = sqrt(0.1)                         # sd
    F_w                       = Normal(μ_w, σ_w)                  # distribution
    min_w::Float64            = -0.5                              # Lower bound
    max_w::Float64            = 1.5                               # Upper bound 
    N_w::Int64                = 41                                # Number of points
    grid_srl_w                = range(min_w, max_w; length = N_w) # Step range
    grid_w::Array{Float64, 1} = collect(grid_srl_w)               # Grid array
    Π_w::Array{Float64, 1}    = diff(append!(cdf(F_w, grid_w), 1))            # probability
end

# Results
mutable struct Results
    W_vf::Array{Float64}       # employed value function
    U_vf::Array{Float64}       # unemployed value function
    S_pf::Array{Float64}       # search policy function
    RW_pf::Array{Float64}      # reservation wage policy function
end

# Initialize results structure
function Initialize()
    @unpack N_h, N_w, T = Primitives()

    # employed workers have 3 state variable 
    # 1 dim is age, 2nd dim is human capital, and 3rd dim is wage
    W_vf = zeros(T, N_h, N_w)

    # unemployed workers have 2 state variable 
    # 1 dim is age, 2nd dim is human capital
    U_vf = zeros(T, N_h)
    S_pf = zeros(T, N_h)
    RW_pf = zeros(T, N_h)

    Results(W_vf, U_vf, S_pf, RW_pf)

end

# search cost function
function c(s::Float64)
    return 0.5*s
end

# probability of drawing an offer
function π(s::Float64)
    return sqrt(s)
end

function vfi!(R::Results; progress::Bool = false)
    @unpack T, β, b, Π_w, grid_w, grid_h, grid_s, ψ_e, ψ_u, δ, N_h = Primitives()

    # terminal period
    R.W_vf[T, :, :] = grid_h * grid_w'
    R.U_vf[T, :] .= b

    # backward induction
    for t=(T-1):-1:1
        if progress
            println(t)
        end

        # employed worker
        for (i_h, h) = enumerate(grid_h), (i_w, w) = enumerate(grid_w)

            # if hc increases, computes the index
            i_h_up = min(i_h+1, N_h)
            
            R.W_vf[t, i_h, i_w] = w * h  # instantaneous utility

            # continuation value
            R.W_vf[t, i_h, i_w] += β * (1 - ψ_e) * (1 - δ) * R.W_vf[t+1, i_h, i_w] # if hc stays same and stays employed
            R.W_vf[t, i_h, i_w] += β * (1 - ψ_e) * δ * R.U_vf[t+1, i_h ]           # if hc stays same and becomes unemployed
            R.W_vf[t, i_h, i_w] += β * ψ_e * (1 - δ) * R.W_vf[t+1, i_h_up, i_w]     # if hc increases and stays employed
            R.W_vf[t, i_h, i_w] += β * ψ_e * δ * R.U_vf[t+1, i_h_up]                # if hc increases and becomes unemployed

        end

        # unemployed worker
        for (i_h, h) = enumerate(grid_h)

            # if hc decreases, computes the index
            i_h_down = max(i_h - 1, 1)

            candidate_max = -Inf

            for (i_s, s) = enumerate(grid_s)

                # instantaneous utility
                candidate_value = b - c(s)

                # continuation value
                # stays unemployed
                candidate_value += β * (1 - ψ_u) * (1 - π(s)) * R.U_vf[t+1, i_h ]  # if hc stays same
                candidate_value += β * ψ_u * (1 - π(s)) * R.U_vf[t+1, i_h_down ]   # if hc decreases

                # gets offer
                candidate_value += β * (1 - ψ_u) * π(s) * Π_w' * max.(R.W_vf[t+1, i_h, :], R.U_vf[t+1, i_h ])    # if hc stays same
                candidate_value += β * ψ_u * π(s) * Π_w' * max.(R.W_vf[t+1, i_h_down, :], R.U_vf[t+1, i_h_down]) # if hc drops

                if candidate_value > candidate_max
                    R.U_vf[t, i_h] = candidate_value
                    R.S_pf[t, i_h] = s
                    R.RW_pf[t, i_h] = minimum(grid_w[R.W_vf[t+1, i_h, :] .> R.U_vf[t+1, i_h]])
                    candidate_max = candidate_value
                end
            end
        end
    end
end