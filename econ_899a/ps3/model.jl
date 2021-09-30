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
using Parameters, CSV, Tables

# Primitive structure
@with_kw struct Primitives
    # life-cycle
    N::Int64                  = 66                               # lifespan
    Jᴿ::Int64                 = 46                               # retirement age
    n::Float64                = 0.011                            # population growth rate

    # preference parameters
    σ::Float64                = 2.0                              # coefficient of relative risk aversion
    β::Float64                = 0.97                             # Discount rate

    # productivity
    η::Array{Float64, 1}      = map(x->parse(Float64,x), readlines("ef.txt")) # deterministic age-efficiency profile
    z_length::Int64           = 2                                # Number of productivity states
    Π::Array{Float64, 2}      = [0.9261 0.0739; 0.0189 0.9811]   # Productivity persistance probabilities
    Π₀::Array{Float64, 1}     = [0.2037, 0.7963]                 # Erodic distribution of Π

    # production
    α::Float64                = 0.36                             # capital elasticity in production
    δ::Float64                = 0.06                             # capital depreciation rate

    # asset grid
    a_min::Float64            = 0.0                              # Lower bound of savings grid
    a_max::Float64            = 75.0                             # Upper bound of savings grid
    a_length::Int64           = 5000                             # Number of points on asset grid
    a_grid_srl                = range(a_min, a_max; length = a_length) # Savings grid step range
    a_grid::Array{Float64, 1} = collect(a_grid_srl)               # Savings grid array
end

mutable struct Results
    θ::Float64                            # proportional labor tax
    z::Array{Float64, 1}                  # Idiosyncratic productivity levels (must be only two states)
    e::Array{Float64, 2}                  # productivity levels
    γ::Float64                            # weight on consumption in utility function
    value_function::Array{Float64}        # value function
    policy_function::Array{Float64}       # savings policy function
    labor_supply::Array{Float64}          # Optimal labor supply
    μ::Array{Float64}                     # Asset distribution
    k::Float64                            # Aggregate capital
    l::Float64                            # Aggregate labor
    w::Float64                            # wage
    r::Float64                            # interest rate
    b::Float64                            # pension benefit
end

function Initialize(θ::Float64, z::Array{Float64, 1}, γ::Float64)
    @unpack a_length, z_length, N, Jᴿ, η = Primitives()

    e = η * z'

    # value_function, policy_function, labor_supply, and F are 3-dimensional arrays
    # 1st dim is age, 2nd dim is asset holding, 3rd dim is productivity state
    value_function       = zeros(N, a_length, z_length)
    policy_function      = zeros(N, a_length, z_length)
    labor_supply         = zeros(Jᴿ-1, a_length, z_length)
    μ                    = ones(N, a_length, z_length) / sum(ones(N, a_length, z_length))


    # Initial values based on representative agent model.
    k = 0.0
    l = 0.0
    w = 1.05
    r = 0.05
    b = 0.2

    Results(θ, z, e, γ, value_function, policy_function, labor_supply, μ, k, l, w, r, b)

end

function Update_prices(results, k, l)
    @unpack α, δ, Jᴿ, N = Primitives()

    results.w = (1 - α) * k ^ α * l ^ (- α)
    results.r = α * k ^ (α - 1) * l ^ (1-α) - δ
    results.b = (results.θ * results.w * l) / sum(results.μ[Jᴿ:N, :, :])
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
    @unpack β, N, Jᴿ, z_length, σ = Primitives()
    @unpack a_min, a_max, a_length, a_grid, a_grid_srl = Primitives()

    # In last period, agents consume everything.
    results.value_function[N,:, 1] = uᴿ.(a_grid, results.γ, σ)

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
                val = uᴿ(budget-a_grid[i_a_p],results.γ,σ)+β*results.value_function[j+1,i_a_p,1]

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
    @unpack β, Π, Jᴿ, z_length, σ = Primitives()
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
                    l = labor_decision(a_grid[i_a], a_grid[i_a_p], results.e[j, i_z], results.θ,
                                           results.γ, results.r, results.w)

                    # Get budget for savings and consumption
                    budget=results.w*(1-results.θ)*results.e[j,i_z]*l + (1+results.r)*a_grid[i_a]

                    val = uᵂ(budget-a_grid[i_a_p], l, results.γ, σ) # Instanteous utility

                    # iterates over tomorrow productivity to add continuation value
                    for i_z_p = 1:z_length
                        val += β*Π[i_z,i_z_p]*results.value_function[j+1,i_a_p,i_z_p]
                    end

                    # if utility starts to decrease, we save values and move on
                    if val < val_previous
                        results.value_function[j, i_a, i_z] = val_previous
                        results.policy_function[j, i_a, i_z] = a_grid[i_a_p-1]
                        results.labor_supply[j,i_a,i_z]=labor_decision(a_grid[i_a],
                                a_grid[i_a_p-1], results.e[j, i_z], results.θ, results.γ, results.r, results.w)
                        choice_lower = i_a_p - 1
                        break

                    # if we're at the top of the grid, we save and move on
                    elseif  i_a_p == a_length
                        results.value_function[j, i_a, i_z] = val
                        results.policy_function[j, i_a, i_z] = a_grid[i_a_p]
                        results.labor_supply[j, i_a, i_z] = labor_decision(a_grid[i_a],
                             a_grid[i_a_p], results.e[j, i_z], results.θ, results.γ, results.r, results.w)
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
################## Functions to solve for asset distribution ###################
################################################################################

function Solve_μ(results::Results; progress::Bool = false)
    @unpack a_grid, N, n, z_length, Π₀, Π, a_length = Primitives()

    # sets distribution to zero.
    results.μ = zeros(N, a_length, z_length)

    # Fills in model-age 1 with erodgic distribution of producitivities.
    results.μ[1, 1, :] = Π₀

    for j = 1:(N-1) # Iterates through model-ages
        if progress
            println(j)
        end
        for i_a in 1:a_length # Iterates through asset levels
            for i_z = 1:z_length
                if results.μ[j, i_a, i_z] == 0 # skips if no mass at j, i_a, i_z
                    continue
                end
                # finds index of assets tomorrow
                i_a_p = argmax(a_grid .== results.policy_function[j, i_a, i_z])
                for i_z_p = 1:z_length # iterates over productivity levels tomorrow
                    results.μ[j+1, i_a_p, i_z_p] += Π[i_z, i_z_p] * results.μ[j, i_a, i_z]
                end
            end
        end
    end

    # sum(results.μ) should equal N at this point (i.e. 1 for each row)
    # Now adjusts for population growth, so that sum(results.μ) = 1

    age_weights_temp = ones(N)

    for i = 1:(N-1)
        age_weights_temp[i + 1] = age_weights_temp[i]/(1+n)
    end

    age_weight = reshape(repeat(age_weights_temp/sum(age_weights_temp), a_length * z_length), N, a_length, z_length)

    results.μ = age_weight .* results.μ
end

################################################################################
######################## Functions for market clearing #########################
################################################################################

function Calculate_labor_supply(results::Results)
    @unpack Jᴿ, a_length, z_length = Primitives()

    e_3d = reshape(repeat(results.e, a_length), Jᴿ -1, a_length, z_length)

    sum(results.μ[1:(Jᴿ - 1),:,:] .* results.labor_supply .* e_3d)
end

function Calculate_capital_supply(results::Results)
    @unpack a_grid, N, z_length, a_length, N = Primitives()

    a_grid_3d = permutedims(reshape(repeat(a_grid, N * z_length), a_length, N, z_length), (2, 1, 3))

    sum(results.μ .* a_grid_3d)
end

function Solve_model(;θ::Float64 = 0.11, z::Array{Float64, 1} = [3.0, 0.5], γ::Float64 = 0.42, λ::Float64 = 0.5)
    results = Initialize(θ, z, γ)

    k_0 = 3.3
    l_0 = 0.3

    Update_prices(results, k_0, l_0)

    ε = 0.001  # tolerence
    i = 0      # counter

    while true
        i += 1
        println("Iteration #", i)

        println("Capital demand: ", k_0)
        println("Labor demand: ", l_0)

        Solve_HH_problem(results)
        Solve_μ(results)

        k_1 = Calculate_capital_supply(results)
        l_1 = Calculate_labor_supply(results)

        println("Capital supply: ", k_1)
        println("Labor supply: ", l_1)

        diff = abs(k_0 - k_1) + abs(l_0 - l_1)

        println("Absolute difference: ", diff)

        println("************************************")

        if diff > ε
            k_0 = λ * k_1 + (1 - λ) * k_0
            l_0 = λ * l_1 + (1 - λ) * l_0
            Update_prices(results, k_0, l_0)
        else
            break
        end
    end

    results.k = k_0
    results.l = l_0

    results
end

################################################################################
############################ Housekeeping function #############################
################################################################################

process_results = function(results::Results)
    # calculate total welfare
    welfare = results.value_function .* results.μ
    welfare = sum(welfare[isfinite.(welfare)])

    # calculate coefficient of variation of wealth
    @unpack a_grid, N, z_length, a_length, N = Primitives()
    a_grid_3d = permutedims(reshape(repeat(a_grid, N * z_length), a_length, N, z_length), (2, 1, 3))
    wealth_mean = sum(results.μ .* a_grid_3d)
    wealth_second_moment = sum(results.μ .* a_grid_3d .^ 2)
    wealth_second_central_moment = wealth_2m - wealth_1m^2
    cv = wealth_mean / sqrt(wealth_second_central_moment)

    # create vector of summary statistics
    [results.θ, results.γ, results.z[1], results.k, results.l,
     results.w, results.r, results.b, welfare, cv]
end

function create_table(results_vector::Array{Results})
    table = DataFrames.DataFrame(Tables.table(reduce(hcat,process_results.(results_vector))'))
    rename!(table, [:theta, :gamma, :z_H, :k, :l, :w, :r, :b, :welfare, :cv])
end
