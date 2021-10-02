# Conesa and Krueger (1999)
# Alex von Hafften
# October 6, 2021

# ECON 899A Computational Economics
# Problem Set 4

# This file contains the functions to estimate transitions between steady states.
# ./run.jl calls this file and runs the model.

include("steady_state.jl")

################################################################################
################## Functions for calculating transition path ###################
################################################################################

mutable struct Transition
    N_t::Int64                           # length of transition
    θ_1::Float64                         # Terminal social security tax level
    k_demand_path::Array{Float64, 1}     # path of capital demand
    l_demand_path::Array{Float64, 1}     # path of labor demand
    w_path::Array{Float64, 1}            # path of wages
    r_path::Array{Float64, 1}            # path of interest rates
    b_path::Array{Float64, 1}            # path of ss benefits
    value_function::Array{Float64}       # 4d array of value functions
    policy_function::Array{Float64}      # 4d array of asset policy functions
    labor_supply::Array{Float64}         # 4d array of optimal labor supply
    μ::Array{Float64}                    # 4d array of asset distribution
end

function Initialize_transition(ss_0::Results, ss_1::Results, N_t::Int64)
    @unpack α, δ, N, a_length, z_length, Jᴿ = Primitives()

    # start with linear transition path
    k_demand_path = collect(ss_0.k .+ ((0:(N_t-1)) ./ (N_t-1)) .* (ss_1.k - ss_0.k))
    l_demand_path = collect(ss_0.l .+ ((0:(N_t-1)) ./ (N_t-1)) .* (ss_1.l - ss_0.l))

    # value_function, policy_function, labor_supply, and μ are 4-dimensional arrays
    # 1st dim is age, 2nd dim is asset holding, 3rd dim is productivity state, 4th dim is transition period
    value_function       = zeros(N, a_length, z_length, N_t)
    policy_function      = zeros(N, a_length, z_length, N_t)
    labor_supply         = zeros(Jᴿ-1, a_length, z_length, N_t)
    μ                    = ones(N, a_length, z_length, N_t) / sum(ones(N, a_length, z_length))

    # Fill in value function, policy function, and labor supply at the end to be the terminal steady state.
    value_function[:, :, :, N_t] = ss_1.value_function
    policy_function[:, :, :, N_t] = ss_1.policy_function
    labor_supply[:, :, :, N_t] = ss_1.labor_supply

    # Fill in μ at the beginning to be that of the initial steady state.
    μ[:, :, :, 1] = ss_1.μ

    # calculate prices
    w_path = (1 - α) .* k_demand_path .^ α .* l_demand_path .^ (- α)
    r_path = α .* k_demand_path .^ (α - 1) .* l_demand_path .^ (1-α) .- δ
    b_path = (ss_1.θ .* w_path .* l_demand_path) ./ reshape(sum(μ[Jᴿ:N, :, :, :], dims = [1, 2, 3]), N_t)

    Transition(N_t, ss_1.θ, k_demand_path, l_demand_path, w_path, r_path, b_path, value_function, policy_function, labor_supply, μ)
end


function Solve_retiree_problem_transition(transition::Transition; progress::Bool = false)
    @unpack β, N, Jᴿ, z_length, σ, γ = Primitives()
    @unpack a_min, a_max, a_length, a_grid, a_grid_srl = Primitives()

    # In last period of life, agents consume everything.
    transition.value_function[N, :, 1, :] = repeat(uᴿ.(a_grid, γ, σ)', transition.N_t)'

    # Backward induction starting at before being in terminal steady state
    for i_t = (transition.N_t-1):-1:1
        if (progress)
            println(i_t)
        end

        # Iterate over ages
        for j = Jᴿ:(N-1)

            choice_lower = 1 # exploits monotonicity of policy function

            # iterates over assets today
            for i_a = 1:a_length

                # budget for savings and consumption
                budget = (1 + transition.r_path[i_t]) * a_grid[i_a] + transition.b_path[i_t]

                val_previous = -Inf # stops grid search when value starts to decrease

                # iterates over assets tomorrow
                for i_a_p = choice_lower:a_length

                    # calculates utility
                    val = uᴿ(budget - a_grid[i_a_p], γ, σ) + β * transition.value_function[j+1,i_a_p,1,i_t+1]

                    # if utility starts to decrease, we save values and move on
                    if val < val_previous
                        transition.value_function[j, i_a, 1, i_t] = val_previous
                        transition.policy_function[j, i_a, 1, i_t] = a_grid[i_a_p-1]
                        choice_lower = i_a_p - 1
                        break

                    # if we're at the top of the grid, we save and move on
                    elseif i_a_p == a_length
                        transition.value_function[j, i_a, 1, i_t] = val
                        transition.policy_function[j, i_a, 1, i_t] = a_grid[i_a_p]
                    end

                    # update val_previous to check if utility is starting to decrease.
                    val_previous = val
                end
            end
        end
    end

    # Fill in for other productivity states (doesn't matter in retirement)
    for i_z = 2:z_length
        transition.policy_function[Jᴿ:N, :, i_z, :] = transition.policy_function[Jᴿ:N, :, 1, :]
        transition.value_function[Jᴿ:N, :, i_z, :] = transition.value_function[Jᴿ:N, :, 1, :]
    end
end


# Solves workers problem. Need to solve retiree problem first.
function Solve_worker_problem_transition(transition::Transition; progress::Bool = false)
    @unpack β, Π, Jᴿ, z_length, σ, γ, e = Primitives()
    @unpack a_min, a_max, a_length, a_grid, a_grid_srl = Primitives()

    # Backward induction starting at terminal steady state.
    for i_t = (transition.N_t-1):-1:1
        if (progress)
            println(i_t)
        end

        # iterates over ages
        for j = 1:(Jᴿ-1)

            # Iterate over productivity today.
            for i_z = 1:z_length

                choice_lower = 1 # exploits monotonicity of policy function

                # Iterate over assets today.
                for i_a = 1:a_length

                    val_previous =  -Inf # stops grid search when value starts to decrease

                    # Iterate over assets tomorrow.
                    for i_a_p = choice_lower:a_length

                        # Solve labor decision (setup only for removal of ss as policy experiment)
                        l = labor_decision(a_grid[i_a], a_grid[i_a_p], e[j, i_z],
                                           transition.θ_1, γ, transition.r_path[i_t],
                                           transition.w_path[i_t])

                        # Get budget for savings and consumption
                        budget = transition.w_path[i_t] * (1.0-transition.θ_1) * e[j,i_z]*l + (1+transition.r_path[i_t])*a_grid[i_a]

                        val = uᵂ(budget-a_grid[i_a_p], l, γ, σ) # Instanteous utility

                        # iterates over tomorrow productivity to add continuation value
                        for i_z_p = 1:z_length
                            val += β * Π[i_z,i_z_p] * transition.value_function[j+1,i_a_p,i_z_p, i_t + 1]
                        end

                        # if utility starts to decrease, we save values and move on
                        if val < val_previous
                            transition.value_function[j, i_a, i_z, i_t] = val_previous
                            transition.policy_function[j, i_a, i_z, i_t] = a_grid[i_a_p-1]
                            transition.labor_supply[j,i_a,i_z, i_t]=labor_decision(a_grid[i_a],
                                       a_grid[i_a_p-1], e[j, i_z], transition.θ_1, γ, transition.r_path[i_t], transition.w_path[i_t])
                            choice_lower = i_a_p - 1
                            break

                        # if we're at the top of the grid, we save and move on
                        elseif  i_a_p == a_length
                            transition.value_function[j, i_a, i_z, i_t] = val
                            transition.policy_function[j, i_a, i_z, i_t] = a_grid[i_a_p]
                            results.labor_supply[j, i_a, i_z, i_t] = labor_decision(a_grid[i_a],
                                a_grid[i_a_p], e[j, i_z], transition.θ_1, γ, transition.r_path[i_t], transition.w_path[i_t])
                        end

                        # update val_previous to check if utility is starting to decrease.
                        val_previous = val
                    end
                end
            end
        end
    end
end


function Solve_HH_problem_transition(transition::Transition)
    Solve_retiree_problem_transition(transition)
    Solve_worker_problem_transition(transition)
end

function Solve_transition(θ_0::Float64, θ_1::Float64,
                          k_0_0::Float64, k_0_1::Float64,
                          l_0_0::Float64, l_0_1::Float64)

    # Steady states
    println("Solve for initial steady state: ")
    ss_0 = Solve_steady_state(k_0_0, l_0_0; θ = θ_0, progress = true)

    println("Solve for terminal steady state: ")
    ss_1 = Solve_steady_state(k_0_1, l_0_1; θ = θ_1, progress = true)

    # Initial guess for transition length
    N_t = 30

    while true # loop for determining length of transition

        transition = Initialize_transition(ss_0, ss_1, N_t)

        while true # loop for convergence of k and l path

            # backward induction to solve HH problem along transition path
            Solve_HH_problem_transition(transition)

            # forward induction to solve μ along transition path

            # calculate k_supply_path and l_supply_path

            # test if k_supply_path and l_supply_path are close to k_demand_path and l_demand_path
            # if not adjust k_demand_path and l_demand_path
        end

        # test if k and l at N are close enough to terminal steady state, ss_1
        # if not increase N

    end

end
