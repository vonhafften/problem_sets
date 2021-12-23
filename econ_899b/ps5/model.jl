# Alex von Hafften 
# Besanko and Doraszelski
# Supplement to Computational Economics
# December 23, 2021

# This file contains the model functions

using Parameters, LinearAlgebra

@with_kw struct Primitives
    # Demand parameters
    a::Float64             = 40.0
    b::Float64             = 10.0

    # Capacity constraint grid
    q_bar_min                  = 0.0
    q_bar_max                  = 45.0
    q_bar_increment            = 5.0
    q_bar_grid::Array{Float64} = q_bar_min:q_bar_increment:q_bar_max
    n_q::Int64                 = length(q_bar_grid)

    # Depreciation (should move to results)
    β::Float64             = 1/1.05
    x_bar::Float64         = 20.0
    θ_bar::Float64         = 0.5
    
end

mutable struct Results

    # Parameters that vary between runs of the model
    competition::String # either quantity or price
    δ::Float64          # depreciation
    α::Float64          # 

    # Static objects
    # Region for bertand
    region::Array{String}

    # optimal quantities at each capacity constraint levels
    q_star_1::Array{Float64}
    q_star_2::Array{Float64}

    # optimal prices at each capacity constraint levels
    p_star_1::Array{Float64}
    p_star_2::Array{Float64}

    # static profit levels at each capacity constraint levels
    π_1::Array{Float64}
    π_2::Array{Float64}

    # Dynamic objects
    # investment policy at each capacity constraint levels
    x_pf_1::Array{Float64}
    x_pf_2::Array{Float64}

    # franchise value at each capacity constraint levels
    vf_1::Array{Float64}
    vf_2::Array{Float64}

    # Stationary distribution
    μ::Array{Float64}
    
end

function Initialize(competition::String, δ::Float64)
    P = Primitives()

    α = P.θ_bar / ((1-δ-P.θ_bar)*P.x_bar)

    region = fill("NA", P.n_q, P.n_q)
    q_star_1 = zeros(P.n_q, P.n_q)
    q_star_2 = zeros(P.n_q, P.n_q)
    p_star_1 = zeros(P.n_q, P.n_q)
    p_star_2 = zeros(P.n_q, P.n_q)
    π_1 = zeros(P.n_q, P.n_q)
    π_2 = zeros(P.n_q, P.n_q)
    x_pf_1 = zeros(P.n_q, P.n_q)
    x_pf_2 = zeros(P.n_q, P.n_q)
    vf_1 = zeros(P.n_q, P.n_q)
    vf_2 = zeros(P.n_q, P.n_q)
    μ = ones(P.n_q, P.n_q) / (P.n_q*P.n_q)

    Results(competition, δ, α, region, q_star_1, q_star_2, p_star_1, p_star_2, π_1, π_2, x_pf_1, x_pf_2, vf_1, vf_2, μ)
end

function demand(price::Float64, P::Primitives)
    return P.a - P.b * price
end

function inv_demand(quantity::Float64, P::Primitives)
    return P.a/P.b - quantity/P.b
end

function compute_static_cournot(q_bar_1::Float64, q_bar_2::Float64, P::Primitives; verbose::Bool = false)
    q_star_1 = q_bar_1
    q_star_2 = q_bar_2
    q_star_1_next = q_bar_1
    q_star_2_next = q_bar_2

    i, maxiter, err = 1, 1000, 100
    i = 1

    while (err > 1e-12) & (maxiter > i)
        # update q_star guesses based on FOC of firm problem with quantity competition
        q_star_1_next = max(0, min(q_bar_1, (P.a - q_star_2) / 2))
        q_star_2_next = max(0, min(q_bar_2, (P.a - q_star_1) / 2))

        err = max(norm(q_star_1_next - q_star_1), norm(q_star_2_next - q_star_2))

        if verbose
            println("Iteration #", i)
            println("Old guess for q_star_1: ", q_star_1)
            println("New guess for q_star_1: ", q_star_1_next)
            println("Old guess for q_star_2: ", q_star_2)
            println("New guess for q_star_2: ", q_star_2_next)
        end

        # update guesses
        q_star_1 = copy(q_star_1_next)
        q_star_2 = copy(q_star_2_next)

        i+=1
        
    end
    return q_star_1, q_star_2
end

function compute_static_bertand(q_bar_1::Float64, q_bar_2::Float64, P::Primitives)
    # optimal quantities in cournot game with zero mc and unlimited capacity
    q_tilde_1 = (P.a - q_bar_2) / 2
    q_tilde_2 = (P.a - q_bar_1) / 2

    # demand at a price of zero
    Q_0 = demand(0.0, P)

    # region A - both firms produce at capacity and sell at the market-clearing price
    if (q_bar_1 <= q_tilde_1) & (q_bar_2 <= q_tilde_2)
        
        price = inv_demand(q_bar_1 + q_bar_2, P)
        return "A", q_bar_1, q_bar_2, price, price, price*q_bar_1, price*q_bar_2

    # region C - both firms each sufficient capacity to serve market.
    # pricing at marginal cost and zero profits.
    elseif (q_bar_1 >= Q_0) & (q_bar_2 >= Q_0) 
        
        # assume that consumer buy equally from firms (it doesn't really matter)
        return "C", Q_0/2, Q_0/2, 0.0, 0.0, 0.0, 0.0
    
    # region B_1 - firm 1 is small than firm 2
    elseif q_bar_1 <= q_bar_2

        # firm 2's problem
        p_star_2 = (P.a - q_bar_1)/(2*P.b)
        q_star_2 = min(q_bar_2, max(0, demand(p_star_2, P) - q_bar_1))
        π_2 = p_star_2 * q_star_2

        # firm 1's problem
        # λ is the multiplier on minimization problem to make better for firm 2 to be the higher priced firm
        λ = sqrt(P.a^2 - 4*P.b*π_2)
        p_2_bar = 1/(2*P.b)*(P.a - λ)

        p_star_1 = p_2_bar
        q_star_1 = min(q_bar_1, demand(p_star_1, P))
        π_1 = p_star_1 * q_star_1
        
        return "B_1", q_star_1, q_star_2, p_star_1, p_star_2, π_1, π_2
    
    # region B_2 - firm 1 is larger than firm 2
    else # q_bar_1 >= q_bar_2
        
        # firm 1's problem
        p_star_1 = (P.a - q_bar_2)/(2*P.b)
        q_star_1 = min(q_bar_1, max(0, demand(p_star_1, P) - q_bar_2))
        π_1 = p_star_1 * q_star_1

        # firm 2's problem
        # λ is the multiplier on minimization problem to make better for firm 1 to be the higher priced firm
        λ = sqrt(P.a^2 - 4*P.b*π_1)
        p_1_bar = 1/(2*P.b)*(P.a - λ)

        p_star_2 = p_1_bar
        q_star_2 = min(q_bar_2, demand(p_star_2, P))
        π_2 = p_star_2 * q_star_2

        return "B_2", q_star_1, q_star_2, p_star_1, p_star_2, π_1, π_2
    
    end
end

function compute_static_results!(R::Results, P::Primitives)

    if R.competition == "quantity"
        for i = 1:P.n_q, j = 1:P.n_q
            R.q_star_1[i,j], R.q_star_2[i,j] = compute_static_cournot(P.q_bar_grid[i], P.q_bar_grid[j], P)
            R.p_star_1[i,j] = inv_demand(R.q_star_1[i,j] + R.q_star_2[i,j], P)
            R.p_star_2[i,j] = inv_demand(R.q_star_1[i,j] + R.q_star_2[i,j], P)
            R.π_1[i,j] = R.p_star_1[i,j] * R.q_star_1[i,j]
            R.π_2[i,j] = R.p_star_2[i,j] * R.q_star_2[i,j]
        end
    elseif R.competition == "price"
        for i = 1:P.n_q, j = 1:P.n_q
            R.region[i,j], R.q_star_1[i,j], R.q_star_2[i,j], R.p_star_1[i,j], R.p_star_2[i,j], R.π_1[i,j], R.π_2[i,j] = compute_static_bertand(P.q_bar_grid[i], P.q_bar_grid[j], P)
        end
    end
end

# returns the probability of going from q_bar to q_bar_next given investment x
function pr(q_bar::Float64, q_bar_next::Float64, x::Float64, δ::Float64, α::Float64, P::Primitives)

    # if at the minimum capacity constraint, you can't drop further down.
    if q_bar == P.q_bar_min
        if q_bar_next < q_bar
            return 0.0
        elseif q_bar_next == q_bar
            return 1/(1+α*x)
        elseif q_bar_next == q_bar + P.q_bar_increment
            return (α*x)/ (1+α*x)
        else
            return 0.0
        end
    # if at the maximum capacity constraint, you can't grow.
    elseif q_bar == P.q_bar_max
        if q_bar_next > q_bar
            return 0.0
        elseif q_bar_next == q_bar
            return (1-δ+α*x)/(1+α*x)/(1+α*x)
        elseif q_bar_next == q_bar - P.q_bar_increment
            return δ/(1+α*x)
        else
            return 0.0
        end
    else # interior capacity constraint levels
        if q_bar_next == q_bar + P.q_bar_increment
            return ((1-δ)*α*x)/(1+α*x)
        elseif q_bar_next == q_bar
            return (1-δ)/(1+α*x) + (δ*α*x)/(1+α*x)
        elseif q_bar_next == q_bar - P.q_bar_increment
            return δ/(1+α*x)
        else
            return 0.0
        end
    end
end

# written from the perspective of firm 1
function W(q_bar_1::Float64, q_bar_2::Float64, x_2::Float64, δ::Float64, α::Float64, vf::Array{Float64}, P::Primitives)
    # get current index of firm 1
    i = Int((q_bar_1 - P.q_bar_min)/P.q_bar_increment + 1)
    if q_bar_1 != P.q_bar_grid[i]
        error("Error in W function")
    end

    result = 0.0

    for j = 1:P.n_q
        result += vf[i,j]*pr(q_bar_2, P.q_bar_grid[j], x_2, δ, α, P)
    end

    return result
end

function Bellman(R::Results, P::Primitives)
    @unpack δ, α = R
    @unpack β, n_q, q_bar_grid, q_bar_increment = P

    # initialize updated guess
    x_pf_1_next = zeros(n_q, n_q)
    x_pf_2_next = zeros(n_q, n_q)
    vf_1_next = zeros(n_q, n_q)
    vf_2_next = zeros(n_q, n_q)

    for i = 1:n_q, j = 1:n_q
        # unpack items
        q_bar_1 = q_bar_grid[i]
        q_bar_2 = q_bar_grid[j]

        # solve firm 1 problem
        x_2 = R.x_pf_2[i, j]
        W_constant = W(q_bar_1,q_bar_2, x_2,δ, α,R.vf_1,P)
        if i == 1
            W_decrease = W_constant
        else
            W_decrease = min(W(q_bar_1-q_bar_increment,q_bar_2, x_2,δ, α,R.vf_1,P), W_constant)
        end
        if i == n_q
            W_increase = W_constant
        else
            W_increase = max(W(q_bar_1+q_bar_increment,q_bar_2, x_2,δ, α,R.vf_1,P), W_constant)
        end

        x_pf_1_next[i,j] = max(0, (-1 + sqrt(β*α*((1-δ)*(W_increase - W_constant) + δ*(W_constant - W_decrease)) ))/α)
        vf_1_next[i,j] = R.π_1[i,j] - x_pf_1_next[i,j]
        for i_p = 1:n_q
            vf_1_next[i,j] += β*W(q_bar_grid[i_p],q_bar_2, x_2,δ, α,R.vf_1,P)*pr(q_bar_1, q_bar_grid[i_p], x_pf_1_next[i,j], δ, α, P)
        end

        # solve firm 2 problem
        vf_2_temp = Float64.(R.vf_2')
        x_1 = R.x_pf_1[i, j]
        W_constant = W(q_bar_2,q_bar_1, x_1,δ, α, vf_2_temp,P)
        if j == 1
            W_decrease = W_constant
        else
            W_decrease = min(W(q_bar_2-q_bar_increment,q_bar_1, x_1,δ, α, vf_2_temp,P), W_constant)
        end
        if j == n_q
            W_increase = W_constant
        else
            W_increase = max(W(q_bar_2+q_bar_increment,q_bar_1, x_1,δ, α, vf_2_temp,P), W_constant)
        end

        x_pf_2_next[i,j] = max(0, (-1 + sqrt(β*α*((1-δ)*(W_increase - W_constant) + δ*(W_constant - W_decrease)) ))/α)
        vf_2_next[i,j] = R.π_2[i,j] - x_pf_2_next[i,j]
        for j_p = 1:n_q
            vf_2_next[i,j] += β * W(q_bar_grid[j_p],q_bar_1, x_1,δ, α, vf_2_temp,P)*pr(q_bar_2, q_bar_grid[j_p], x_pf_2_next[i,j], δ, α, P)
        end     
    end
    x_pf_1_next, x_pf_2_next, vf_1_next, vf_2_next
end

function Solve_model(competition::String, δ::Float64; verbose::Bool = false)
    R = Initialize(competition, δ)
    P = Primitives()
    
    compute_static_results!(R, P)

    i, maxiter, err = 1, 1000, 100
    i = 1
    
    while (err > 1e-12) & (maxiter > i)
        x_pf_1_next, x_pf_2_next, vf_1_next, vf_2_next = Bellman(R, P)

        err = max(norm(x_pf_1_next - R.x_pf_1), norm(x_pf_2_next - R.x_pf_2), norm(vf_1_next - R.vf_1), norm(vf_2_next - R.vf_2))

        if verbose
            println("Iteration #", i)
            println("Sup norm: ", err)
        end

        R.x_pf_1 = copy(x_pf_1_next)
        R.x_pf_2 = copy(x_pf_2_next)
        R.vf_1 = copy(vf_1_next)
        R.vf_2 = copy(vf_2_next)

        i+=1
    end

    println("Converged in ", i, " iterations.")
    
    return R
end

