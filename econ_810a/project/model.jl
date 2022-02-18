# Alex von Hafften
# ECON 810A - Project
# Bewley with long-term illiquid assets

using Parameters

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_810a/project/")

# Primitive structure
@with_kw struct Primitives

    # parameters
    δ::Float64                = 0.8                               # stochastic maturity
    r_S::Float64              = 0.01                              # exogeneous short-term interest rate
    r_L::Float64              = 0.09                              # exogeneous long-term interest rate
    κ::Float64                = 0.5                               # firesale penalty rate
    β::Float64                = (1/(1+r_S))^(1/12)                # Discount rate
    γ::Float64                = 2.0                               # Coefficient of relative risk aversion

    # asset grid
    min_a::Float64            = 0.0                               # Lower bound
    max_a::Float64            = 5.0                              # Upper bound
    N_a::Int64                = 100                               # Number of points
    grid_srl_a                = range(min_a, max_a; length = N_a) # Step range
    grid_a::Array{Float64, 1} = collect(grid_srl_a)               # Grid array

    # liquidiate grid
    min_l::Float64            = 0.0                               # Lower bound
    max_l::Float64            = 1.0                               # Upper bound
    N_l::Int64                = 101                                # Number of points
    grid_srl_l                = range(min_l, max_l; length = N_l) # Step range
    grid_l::Array{Float64, 1} = collect(grid_srl_l)               # Grid array

    # exogenous labor earnings process
    grid_y::Array{Float64, 1} = [2.0, 0.5]           # exogeneous labor earning states
    Π_y::Array{Float64, 2}    = [0.9 0.1; 0.5 0.5] # transition probabilities
    N_y::Int64                = 2                    # Number of points
end

# Results
mutable struct Results
    vf::Array{Float64}         # value function
    pf_c::Array{Float64}       # consumption policy function
    pf_a::Array{Float64}       # asset policy function
    pf_l::Array{Float64}       # liquidation policy function
end

# Initialize results structure
function Initialize()
    @unpack N_y, N_a = Primitives()

    # HH has 2 state variable 
    # 1 dim is assets, 2nd dim is employment status
    vf   = zeros(N_a, N_y)
    pf_c = zeros(N_a, N_y)
    pf_a = zeros(N_a, N_y)
    pf_l = zeros(N_a, N_y)

    Results(vf, pf_c, pf_a, pf_l)
end

# utility function
function u(c::Float64, γ::Float64)
    if (c > 0)
        return (c^(1-γ))/(1-γ)
    else
        return -1/eps()
    end
end

function fire_sale(a::Float64)
    return log(1 + a)
end

################################################################################
######################### Functions to solve HH problem ########################
################################################################################

function HH_Bellman(R::Results; progress::Bool = false)
    @unpack β, γ, Π_y, grid_y, grid_a, grid_l, δ, κ, r_S, r_L, N_a, N_y = Primitives()

    vf_next = zeros(N_a, N_y)  # Initialize next value function iteration

    # cycle over employment statuses and asset today
    for (i_y, y) = enumerate(grid_y), (i_a, a) = enumerate(grid_a)

        if progress
            println(y, ", ", a)
        end
        
        candidate_max = -1/eps()

        # Iterate over fractions of your long-term assets to liquidate l
        for (i_l, l) = enumerate(grid_l)

            # budget equals labor earnings + maturing assets + fire-sold long-term assets
            budget = y + (1+r_S)*δ*a + fire_sale(l*(1+r_L)*(1-δ)*a)
            lt_assets = (1-l)*(1+r_L)*(1-δ)*a

            # if you consume nothing the max asset position next period
            max_a_p = budget + lt_assets

            # if you consume your entire budget
            min_a_p = lt_assets

            # limit asset grid to be feasible choices
            grid_a_p = grid_a[(grid_a .> min_a_p).*(grid_a .< max_a_p)]

            for (i_a_p, a_p) = enumerate(grid_a_p)
                consumption = budget - (a_p - lt_assets)

                value = u(consumption, γ) + β * Π_y[i_y, :]' * R.vf[i_a_p, :]

                if value > candidate_max
                    vf_next[i_a, i_y] = value
                    R.pf_c[i_a, i_y] = consumption
                    R.pf_a[i_a, i_y] = a_p
                    R.pf_l[i_a, i_y] = l
                    candidate_max = value
                end
            end
        end
    end
    vf_next
end


function Solve_HH_problem!(R::Results, tolerence::Float64 = 1e-5; progress::Bool = false)

    err, i, max_iter = 100.0, 0, 100

    while err > tolerence
        i += 1

        vf_next = HH_Bellman(R)
        err = maximum(abs.(vf_next.-R.vf))
        R.vf = vf_next

        if progress
            println("Iteration: ", i)
            println("Error: ", err)
        end

        if i > max_iter
            println("Maximum iterations reached.")
            break
        end
    end

    println("HH Problem converged in ", i, " iterations")
end


using Plots

R = Initialize()
Solve_HH_problem!(R; progress = true)

# plots 

@unpack grid_a, δ = Primitives()

plot(grid_a, R.vf, label = ["Employed" "Unemployed"], legend = :bottomright);
title!("Value Function");
xlabel!("Assets")
savefig("vf.png")

plot(grid_a, R.pf_c, label = ["Employed" "Unemployed"], legend = :bottomright);
title!("Consumption Policy Function");
xlabel!("Assets")
savefig("consumption_pf.png")

plot(grid_a, R.pf_a, label =["Employed" "Unemployed"], legend = :bottomright);
plot!(grid_a, grid_a, label = "45 degree line");
title!("Asset Policy Function");
xlabel!("Assets")
savefig("asset_pf.png")

plot(grid_a, R.pf_l, label =["Employed" "Unemployed"], legend = :bottomright);
title!("Liquiation Policy Function");
xlabel!("Assets")
savefig("liquidation_pf.png")

plot(grid_a, (1 .-R.pf_l).*[grid_a grid_a].*(1-δ), label =["Employed" "Unemployed"], legend = :bottomright);
plot!(grid_a, grid_a*(1-δ), label = "Total Long-Term Assets");
title!("Long-Term Assets Function");
xlabel!("Assets")
savefig("liquidated_asset_pf.png")


plot(grid_a[:,1], R.pf_c[:,1], legend = :bottomright, label = "consumption")
plot!(grid_a[:,1], R.pf_c[:,1] .+ R.pf_a[:,1], legend = :bottomright, label = "consumption + savings")
plot!(grid_a, grid_a .+ 2, label = "Total")
plot!(grid_a, grid_a .+ 2, label = "45 degree line")
title!("Consumption Policy Function");
xlabel!("Assets")