using Parameters, Interpolations, Plots

@with_kw struct Primitives
    α::Float64              = 1.5                                        # coefficient of relative risk aversion
    β::Float64              = 0.8                                        # discount factor
    r::Float64              = 0.04                                       # real risk-free interest rate
    ρ::Float64              = 0.9                                        # legal record keeping technology parameter
end

@with_kw struct Grids
    # earning states
    grid_y::Vector{Float64} = [0.05, 1.0]                                # income states
    N_y::Int64              = length(grid_y)                             # number of earning states
    Π_y::Matrix{Float64}    = [0.75 0.25; 0.25 0.75]                     # transition probabilities

    # assets
    min_a::Float64          = -1.0                                       # minimum a
    max_a::Float64          = 4.0                                        # maximum a
    N_a::Int64              = 501                                        # number of a points
    grid_a::Vector{Float64} = collect(range(min_a, max_a; length = N_a)) # a grid

    # saving
    grid_a_p::Vector{Float64} = grid_a[grid_a .>= 0.0]                   # a grid
    N_a_p::Int64            = length(grid_a_p)                           # number of a points

    # borrowing
    grid_a_n::Vector{Float64} = grid_a[grid_a .< 0.0]                    # a grid
    N_a_n::Int64            = length(grid_a_n)                           # number of a points
end

mutable struct Results
    vf_h0::Array{Float64}   # value function for agents without flag
    pf_a_h0::Array{Float64} # policy function for saving/borrowing for agents without flag
    pf_d_h0::Array{Float64} # policy function for defaulting for agents without flag
    vf_h1::Array{Float64}   # value function for agents with flag
    pf_a_h1::Array{Float64} # policy function for saving for agents with flag
    q::Matrix{Float64}      # discount bond prices
end

# Initialize Results struct
function Initialize_Results(P::Primitives, G::Grids)
    vf_h0   = zeros(G.N_y, G.N_a)
    pf_a_h0 = zeros(G.N_y, G.N_a)
    pf_d_h0 = fill(0, G.N_y, G.N_a)
    vf_h1   = zeros(G.N_y, G.N_a_p)
    pf_a_h1 = zeros(G.N_y, G.N_a_p)
    q       = zeros(G.N_y, G.N_a_n) .+ 1/(1+P.r)

    Results(vf_h0, pf_a_h0, pf_d_h0, vf_h1, pf_a_h1, q)
end

# flow utility function
function u(c::Float64, P::Primitives)
    if c >= 0.0
        return (c^(1-P.α))/(1-P.α)
    else
        return -Inf
    end
end

# bellman operator for agent without bankruptcy flag
function apply_bellman_h0!(P::Primitives, G::Grids, R::Results)
    vf_h0_next = zeros(G.N_y, G.N_a)

    for (i_y, y) = enumerate(G.grid_y),  (i_a, a) = enumerate(G.grid_a)

        # how would default look?
        vf_h0_next_d_candidate = u(y, P)
        for (i_y_p, y_p) = enumerate(G.grid_y)
            vf_h0_next_d_candidate += P.β * G.Π_y[i_y, i_y_p] * R.vf_h1[i_y_p, 1]
        end

        # how does not defaulting look?
        vf_h0_next_nd_candidate = -Inf
        for (i_a_p, a_p) = enumerate(G.grid_a)
            if a_p < 0.0
                value = u(y + a - R.q[i_y, i_a_p] * a_p, P)
            else
                value = u(y + a - 1/(1+P.r) * a_p, P)
            end
            for (i_y_p, y_p) = enumerate(G.grid_y)
                value += P.β * G.Π_y[i_y, i_y_p] * R.vf_h0[i_y_p, i_a_p]
            end
            if value > vf_h0_next_nd_candidate
                vf_h0_next_nd_candidate = value
                R.pf_a_h0[i_y, i_a] = a_p
            end
        end

        # better to default or not default?
        if vf_h0_next_d_candidate > vf_h0_next_nd_candidate
            vf_h0_next[i_y, i_a] = vf_h0_next_d_candidate
            R.pf_a_h0[i_y, i_a] = minimum(G.grid_a_p)
            R.pf_d_h0[i_y, i_a] = 1
        else
            vf_h0_next[i_y, i_a] = vf_h0_next_nd_candidate
            R.pf_d_h0[i_y, i_a] = 0
        end
    end

    vf_h0_next
end

function apply_bellman_h1!(P::Primitives, G::Grids, R::Results)
    
    vf_h1_next = zeros(G.N_y, G.N_a_p)

    # with bankruptcy flag
    for (i_y, y) = enumerate(G.grid_y),  (i_a, a) = enumerate(G.grid_a_p)
        vf_h1_next_candidate = -Inf
        for (i_a_p, a_p) = enumerate(G.grid_a_p)
            value = u(y + a - 1/(1+P.r) * a_p, P)
            for (i_y_p, y_p) = enumerate(G.grid_y)
                value += P.β * G.Π_y[i_y, i_y_p] * P.ρ * R.vf_h1[i_y_p, i_a_p]
                value += P.β * G.Π_y[i_y, i_y_p] * (1 - P.ρ) * R.vf_h0[i_y_p, i_a_p + G.N_a - G.N_a_p ]
            end
            if value > vf_h1_next_candidate
                vf_h1_next[i_y, i_a] = value
                vf_h1_next_candidate = value
                R.pf_a_h1[i_y, i_a] = a_p
            end
        end
    end
    vf_h1_next
end

function solve_vf!(P::Primitives, G::Grids, R::Results)
    tol = 0.001
    err = 100.0
    max_iteration = 1000
    i = 1

    while true
        vf_h0_next = apply_bellman_h0!(P, G, R)
        vf_h1_next = apply_bellman_h1!(P, G, R)
        err = max(maximum(abs.(R.vf_h0 - vf_h0_next)), maximum(abs.(R.vf_h1 - vf_h1_next)))
        R.vf_h0 = vf_h0_next
        R.vf_h1 = vf_h1_next
        # println(err)
        if err < tol
            break
        end
        if i > max_iteration
            break
        end
        i += 1
    end
end

function compute_q(P::Primitives, G::Grids, R::Results)
    q_next = zeros(G.N_y, G.N_a_n) .+ 1/(1+P.r)
    
    for (i_y, y) = enumerate(G.grid_y),  (i_a_p, a_p) = enumerate(G.grid_a_n)
        delta = 0
        for (i_y_p, y_p) = enumerate(G.grid_y)
            delta += G.Π_y[i_y, i_y_p] * R.pf_d_h0[i_y_p, i_a_p]
        end 
        q_next[i_y, i_a_p] = (1 - delta) / (1+P.r)
    end
    q_next
end

function solve_q!(P::Primitives, G::Grids, R::Results)

    tol = 0.001
    err = 100.0
    max_iteration = 1000
    i = 1

    λ = 0.5

    while true
        display(plot(G.grid_a_n, R.q', label = "y = " .* string.(G.grid_y'), title = "q", color = ["red" "blue"], legend = :bottomright))

        solve_vf!(P, G, R)
        q_next = compute_q(P, G, R)

        err = maximum(abs.(R.q - q_next))
        R.q = λ*R.q + (1-λ)*q_next

        println(err)
        if err < tol
            break
        end
        if i > max_iteration
            break
        end
        i += 1
    end
end
