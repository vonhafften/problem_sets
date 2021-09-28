# Conesa and Krueger (1999)
# Alex von Hafften
# September 22, 2021

# ECON 899A Computational Economics
# Problem set 3

# This file contains the functions to estimate the model.
# ./run.jl calls this file and runs the model.

################################################################################
######################### Housekeeping functions ###############################
################################################################################

# Load libraries
using Parameters

# Primitive structure
@with_kw struct Primitives
    N::Int64                  = 66                               # lifespan
    n::Float64                = 0.011                            # population growth rate
    Jᴿ::Int64                 = 46                               # retirement age
    θ::Float64                = 0.11                             # proportional labor tax
    γ::Float64                = 0.42                             # weight on consumption in utility function
    σ::Float64                = 2.0                              # coefficient of relative risk aversion
    η::Array{Float64, 1}      = map(x->parse(Float64,x), readlines("ef.txt")) # deterministic age-efficiency profile
    z::Array{Float64, 1}      = [3.0, 0.5]                       # Idiosyncratic productivity levels
    z_length::Int64           = length(z)                        # Number of productivity states
    Π::Array{Float64, 2}      = [0.9261 0.0739; 0.0189 0.9811]   # Productivity persistance probabilities
    Π₀::Array{Float64, 1}     = [0.2037, 0.7963]                 # Erodic distribution of Π
    e::Array{Float64, 2}      = η * z'                           # productivity levels
    α::Float64                = 0.36                             # capital elasticity in production
    δ::Float64                = 0.06                             # capital depreciation rate
    β::Float64                = 0.97                             # Discount rate
    a_min::Float64            = 0.0                              # Lower bound of savings grid
    a_max::Float64            = 75.0                             # Upper bound of savings grid
    a_length::Int64           = 5000                             # Number of points on asset grid
    a_grid_srl                = range(a_min, a_max; length = a_length) # Savings grid step range
    a_grid::Array{Float64, 1} = collect(a_grid_srl)               # Savings grid array
end

mutable struct Results
    value_function::Array{Float64}        # value function
    policy_function::Array{Float64}       # savings policy function
    labor_supply::Array{Float64}          # Optimal labor supply
    F::Array{Float64}                     # Asset distribution
    w::Float64                            # wage
    r::Float64                            # interest rate
    b::Float64                            # pension benefit
end

function Initialize()
    @unpack a_length, z_length, N, Jᴿ = Primitives()

    # value_function, policy_function, labor_supply, and F are 3-dimensional arrays
    # 1st dim is age, 2nd dim is asset holding, 3rd dim is productivity state
    value_function       = zeros(N, a_length, z_length)
    policy_function      = zeros(N, a_length, z_length)
    labor_supply         = zeros(Jᴿ-1, a_length, z_length)
    F                    = zeros(N, a_length, z_length)

    # Initial values based on representative agent model.
    w = 1.05
    r = 0.05
    b = 0.2

    Results(value_function, policy_function, labor_supply, F, w, r, b)
end

################################################################################
######################### Functions to solve HH problem ########################
################################################################################

# Utility function of retired agent
function uᴿ(c::Float64, γ::Float64, σ::Float64)
    if (c > 0)
        (c^((1-σ) * γ))/(1 - σ)
    else
        -Inf
    end
end

function Solve_retiree_problem(results::Results; progress::Bool = false)
    @unpack β, N, Jᴿ, z_length, e, γ, σ = Primitives()
    @unpack a_min, a_max, a_length, a_grid, a_grid_srl = Primitives()

    # In last period, agents consume everything.
    results.value_function[N,:, 1] = uᴿ.(a_grid, γ, σ)

    # Backward induction starting at period before last period
    for j = (N-1):-1:Jᴿ
        if (progress)
            println(j)
        end

        choice_lower = 1 # exploits monotonicity of policy function

        # iterates over assets today
        for i_a = 1:a_length

            # budget for savings and consumption
            budget = (1 + results.r) * a_grid[i_a] + results.b

            val_previous = -Inf # stops grid search when value starts to decrease

            # iterates over assets tomorrow
            for i_a_p = choice_lower:a_length

                # calculates utility
                val = uᴿ(budget-a_grid[i_a_p],γ,σ)+β*results.value_function[j+1,i_a_p,1]

                # if utility starts to decrease, we save values and move on
                if val < val_previous
                    results.value_function[j, i_a, 1] = val_previous
                    results.policy_function[j, i_a, 1] = a_grid[i_a_p-1]
                    choice_lower = i_a_p - 1
                    break

                # if we're at the top of the grid, we save and move on
                elseif i_a_p == a_length
                    results.value_function[j, i_a, 1] = val
                    results.policy_function[j, i_a, 1] = a_grid[i_a_p]
                end

                # update val_previous to check if utility is starting to decrease.
                val_previous = val
            end
        end
    end

    # Fill in for other productivity states (doesn't matter in retirement)
    for i_z = 2:z_length
        results.policy_function[:, :, i_z] = results.policy_function[:, :, 1]
        results.value_function[:, :, i_z] = results.value_function[:, :, 1]
    end
end

# static labor supply decision
function labor_decision(a::Float64, a_p::Float64, e::Float64, θ::Float64,
                        γ::Float64, r::Float64, w::Float64)
    interior_solution = (γ*(1-θ)*e*w - (1-γ)*((1+r)*a - a_p)) / ((1-θ)*w*e)
    min(1, max(0, interior_solution))
end

# utility function of workers
function uᵂ(c::Float64, l::Float64, γ::Float64, σ::Float64)
    if (c > 0 && l >= 0 && l <= 1)
        (((c^γ) * ((1 - l)^(1-γ)))^(1-σ))/(1 - σ)
    else
        -Inf
    end
end

# Solves workers problem. Need to solve retiree problem first.
function Solve_worker_problem(results::Results; progress::Bool = false)
    @unpack β, Π, Jᴿ, θ, e, z_length, γ, σ = Primitives()
    @unpack a_min, a_max, a_length, a_grid, a_grid_srl = Primitives()

    # Backward induction starting period before retirement.
    for j = (Jᴿ-1):-1:1
        if (progress)
            println(j)
        end

        # Iterate over productivity today.
        for i_z = 1:z_length

            choice_lower = 1 # exploits monotonicity of policy function

            # Iterate over assets today.
            for i_a = 1:a_length

                val_previous = -Inf # stops grid search when value starts to decrease

                # Iterate over assets tomorrow.
                for i_a_p = choice_lower:a_length

                    # Solve labor decision
                    l = labor_decision(a_grid[i_a], a_grid[i_a_p], e[j, i_z], θ,
                                           γ, results.r, results.w)

                    # Get budget for savings and consumption
                    budget=results.w*(1-θ)*e[j,i_z]*l + (1+results.r)*a_grid[i_a]

                    val = uᵂ(budget-a_grid[i_a_p], l, γ, σ) # Instanteous utility

                    # iterates over tomorrow productivity to add continuation value
                    for i_z_p = 1:z_length
                        val += β*Π[i_z,i_z_p]*results.value_function[j+1,i_a_p,i_z_p]
                    end

                    # if utility starts to decrease, we save values and move on
                    if val < val_previous
                        results.value_function[j, i_a, i_z] = val_previous
                        results.policy_function[j, i_a, i_z] = a_grid[i_a_p-1]
                        results.labor_supply[j,i_a,i_z]=labor_decision(a_grid[i_a],
                                a_grid[i_a_p-1], e[j, i_z], θ, γ, results.r, results.w)
                        choice_lower = i_a_p - 1
                        break

                    # if we're at the top of the grid, we save and move on
                    elseif  i_a_p == a_length
                        results.value_function[j, i_a, i_z] = val
                        results.policy_function[j, i_a, i_z] = a_grid[i_a_p]
                        results.labor_supply[j, i_a, i_z] = labor_decision(a_grid[i_a],
                             a_grid[i_a_p], e[j, i_z], θ, γ, results.r, results.w)
                    end

                   # update val_previous to check if utility is starting to decrease.
                    val_previous = val
                end
            end
        end
    end
end

# Solves HH problem. Retiree first, then worker.
function Solve_HH_problem(results::Results; progress::Bool = false)

    Solve_retiree_problem(results)
    Solve_worker_problem(results)

end

################################################################################
######## Functions to solve for invariant asset distribution ###################
################################################################################

function Solve_F(results::Results)
    @unpack N, n, z_length, Π₀, a_length = Primitives()

    μ = ones(N)

    for i = 1:(N-1)
        μ[i + 1] = μ[i]/(1+n)
    end

    μ = μ/sum(μ)

    for i = 1:z_length
        results.F[1, 1, i] = μ[1] * Π₀[i]
    end

    for j = 2:N
        for i_a = 1:a_length
            for i_z = 1:z_length
                results.policy_function[j, i_a, i_z]

                ###############################################################
                ###############################################################
                results.F[1, 1, i]
            end
        end
    end
end
