# Alex von Hafften
# ECON 810A - Project
# Bewley with long-term illiquid assets
# this code was my first attempt. Liquidate was a binary choice. 
# Not very interesting results and it took forever to converge. Moving to continuous liquidation choice.

using Parameters

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_810a/project/")

# Primitive structure
@with_kw struct Primitives

    # parameters
    δ::Float64                = 0.2                               # stochastic maturity
    r_S::Float64              = 0.03                              # exogeneous short-term interest rate
    r_L::Float64              = 0.05                              # exogeneous long-term interest rate
    κ::Float64                = 0.8                               # discount rate
    β::Float64                = (1/(1+r_S))^(1/12)               # Discount rate
    γ::Float64                = 2.0                               # Coefficient of relative risk aversion

    # asset grid
    min_a::Float64            = 0.00                              # Lower bound
    max_a::Float64            = 10.0                              # Upper bound
    N_a::Int64                = 100                               # Number of points
    grid_srl_a                = range(min_a, max_a; length = N_a) # Step range
    grid_a::Array{Float64, 1} = collect(grid_srl_a)               # Grid array

    # exogenous labor earnings process
    grid_y::Array{Float64, 1} = [1, 0.5]             # exogeneous labor earning states
    Π_y::Array{Float64, 2}    = [0.97 0.03; 0.5 0.5] # transition probabilities
    N_y::Int64                = 2                    # Number of points
end

# Results
mutable struct Results
    vf::Array{Float64}         # value function
    pf_a::Array{Float64}       # asset policy function
    pf_l::Array{Int64}        # liquidation policy function
end

# Initialize results structure
function Initialize()
    @unpack N_y, N_a = Primitives()

    # HH has 2 state variable 
    # 1 dim is assets, 2nd dim is employment status
    vf = zeros(N_a, N_y)
    pf_a = zeros(N_a, N_y)
    pf_l = fill(0, N_a, N_y)

    Results(vf, pf_a, pf_l)
end

# utility function
function u(c::Float64, γ::Float64)
    if (c > 0)
        return (c^(1-γ))/(1-γ)
    else
        return -1/eps()
    end
end


################################################################################
######################### Functions to solve HH problem ########################
################################################################################

function HH_Bellman(R::Results; progress::Bool = false)
    @unpack β, γ, Π_y, grid_y, grid_a, δ, κ, r_S, r_L, N_a, N_y = Primitives()

    vf_next = zeros(N_a, N_y)  # Initialize next value function iteration

    # cycle over employment statuses and asset today
    for (i_y, y) = enumerate(grid_y), (i_a, a) = enumerate(grid_a)

        if progress
            println(y, ", ", a)
        end

        # if not liquidiating long-term assets
        candidate_max_nliq = -1/eps()
        candidate_a_p_nliq = 0.0
        budget = y + δ*a*(1+r_S)
        for (i_a_d_p, a_d_p) = enumerate(grid_a)
            a_p = (1-δ)*a*(1+r_L) + a_d_p
            i_a_p = argmin(abs.(grid_a .- a_p))
            value = u(budget - a_d_p, γ) + β * Π_y[i_y, :]' * R.vf[i_a_p, :]
            if value > candidate_max_nliq
                candidate_max_nliq = value
                candidate_a_p_nliq = a_p
            end
        end

        # if liquidiating long-term assets
        candidate_max_liq = -1/eps()
        candidate_a_p_liq = 0.0
        budget = y + δ*a*(1+r_S) + κ*(1-δ)*a*(1+r_L)
        for (i_a_p, a_p) = enumerate(grid_a)
            value = u(budget - a_p, γ) + β * Π_y[i_y, :]' * R.vf[i_a_p, :]
            if value > candidate_max_liq
                candidate_max_liq = value
                candidate_a_p_liq = a_p
            end
        end

        if candidate_max_nliq >= candidate_max_liq
            vf_next[i_a, i_y] = candidate_max_nliq
            R.pf_a[i_a, i_y] = candidate_a_p_nliq
            R.pf_l[i_a, i_y] = 0
        elseif candidate_max_nliq < candidate_max_liq
            vf_next[i_a, i_y] = candidate_max_liq
            R.pf_a[i_a, i_y] = candidate_a_p_liq
            R.pf_l[i_a, i_y] = 1
        else
            error("Issue with Bellman")
        end
    end
    vf_next
end


function Solve_HH_problem!(R::Results, tolerence::Float64 = 1e-5; progress::Bool = false)

    err, i, max_iter = 100.0, 0, 1000

    while err > tolerence
        i += 1

        if progress
            println("Iteration: ", i)
            println("Error: ", err)
        end

        vf_next = HH_Bellman(R)
        err = abs.(maximum(vf_next.-R.vf))
        R.vf = vf_next

        if i > max_iter
            println("Maximum iterations reached.")
            break
        end
    end

    println("HH Problem converged in ", i, " iterations")
end

R = Initialize()
Solve_HH_problem!(R;progress = true)
plot(R.vf)
plot(R.pf_a)
plot(R.pf_l)