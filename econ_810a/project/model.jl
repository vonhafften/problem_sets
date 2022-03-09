# Alex von Hafften
# ECON 810A - Project
# Bewley with long-term illiquid assets

using Parameters, Interpolations, Optim, Roots

cd("/Users/vonhafften/Documents/UW Madison/problem_sets/econ_810a/project/")

# Primitive structure
@with_kw struct Primitives

    # parameters
    δ::Float64                = 0.17143                           # stochastic maturity
    r::Float64                = 0.1                               # exogeneous LT bond interest rate
    β::Float64                = 0.99                              # discount rate
    σ::Float64                = 2.0                               # coefficient of relative risk aversion
    γ::Float64                = 0.2                               # probability of meeting intermediary
    T::Int64                  = 30                                # lifespan
    θ::Float64                = 0.5                               # hh bargaining power

    # cash grid
    min_a::Float64            = 0.0                               # Lower bound
    max_a::Float64            = 8.0                               # Upper bound
    N_a::Int64                = 21                                # Number of points
    grid_srl_a                = range(min_a, max_a; length = N_a) # Step range
    grid_a::Array{Float64, 1} = collect(grid_srl_a)               # Grid array

    # LT bond grid
    min_b::Float64            = 0.0                               # Lower bound
    max_b::Float64            = 8.0                               # Upper bound
    N_b::Int64                = 21                                # Number of points
    grid_srl_b                = range(min_b, max_b; length = N_b) # Step range
    grid_b::Array{Float64, 1} = collect(grid_srl_b)               # Grid array

    # liquidiate grid
    min_l::Float64            = 0.1                               # Lower bound
    max_l::Float64            = 1.0                               # Upper bound
    N_l::Int64                = 10                                # Number of points
    grid_srl_l                = range(min_l, max_l; length = N_l) # Step range
    grid_l::Array{Float64, 1} = collect(grid_srl_l)               # Grid array

    # exogenous labor earnings process
    grid_y::Array{Float64, 1} = [1.0, 0.5]           # exogeneous labor earning states
    Π_y::Array{Float64, 2}    = [0.9 0.1; 0.5 0.5]   # transition probabilities
    N_y::Int64                = 2                    # Number of points
end

# Results
mutable struct Results

    β_I::Float64               # discount rate of intermediary
    γ::Float64                 # probability of meeting intermediary

    vf_m::Array{Float64}       # matched value function
    vf_u::Array{Float64}       # unmatched value function

    pf_m_c::Array{Float64}     # matched consumption policy function
    pf_u_c::Array{Float64}     # unmatched consumption policy function
    
    pf_m_a::Array{Float64}     # matched cash policy function
    pf_u_a::Array{Float64}     # unmatched cash policy function

    pf_m_b_tilde::Array{Float64} # matched LT bond purchase policy function
    pf_u_b_tilde::Array{Float64} # unmatched LT bond purchase policy function

    pf_m_b_hat::Array{Float64}   # matched LT bond sales policy function
    pf_m_l::Array{Float64}       # matched liquidation policy function
    price::Array{Float64}        # price of LT bond sale
end

# Initialize results structure
function Initialize(β_I::Float64, γ::Float64)
    @unpack T, N_y, N_a, N_b = Primitives()

    # HH has 2 state variable 
    # 1 dim is age, 2 dim is employment status, 3 dim is cash, 4 dim is LT bonds
    vf_m         = zeros(T, N_y, N_a, N_b)
    vf_u         = zeros(T, N_y, N_a, N_b)
    pf_m_c       = zeros(T, N_y, N_a, N_b)
    pf_u_c       = zeros(T, N_y, N_a, N_b)
    pf_m_a       = zeros(T, N_y, N_a, N_b)
    pf_u_a       = zeros(T, N_y, N_a, N_b)
    pf_m_b_tilde = zeros(T, N_y, N_a, N_b)
    pf_u_b_tilde = zeros(T, N_y, N_a, N_b)
    pf_m_b_hat   = zeros(T, N_y, N_a, N_b)
    pf_m_l       = zeros(T, N_y, N_a, N_b)
    price        = zeros(T, N_y, N_a, N_b)

    Results(β_I, γ, vf_m, vf_u, pf_m_c, pf_u_c, pf_m_a, pf_u_a, pf_m_b_tilde, pf_u_b_tilde, pf_m_b_hat, pf_m_l, price)
end

# utility function
function u(c::Float64, σ::Float64)
    if (c > 0)
        return (c^(1-σ))/(1-σ)
    else
        return -1/eps()
    end
end

function W(b::Float64, δ::Float64, β_I::Float64)
    return ((β_I * δ) / (1- β_I*(1-δ))) * b
end

################################################################################
######################### Functions to solve HH problem ########################
################################################################################

function Solve_terminal_period!(R::Results)
    P = Primitives()

    for (i_y, y) = enumerate(P.grid_y), (i_a, a) = enumerate(P.grid_a), (i_b, b) = enumerate(P.grid_b)
        
        # all agents consume all their labor earnings, cash, and matured LT bonds
        consumption = y + a + P.δ * b * (1 + P.r)
        R.vf_u[P.T, i_y, i_a, i_b] = u(consumption, P.σ)
        R.vf_m[P.T, i_y, i_a, i_b] = u(consumption, P.σ)

    end    
    println(P.T)
end

function Solve_nonterminal_period!(R::Results)
    P = Primitives()

    # iterate backward
    for t = (P.T-1):-1:1

        println(t)

        vf_m_e_i = LinearInterpolation((P.grid_a, P.grid_b), R.vf_m[t+1, 1, :, :])
        vf_m_u_i = LinearInterpolation((P.grid_a, P.grid_b), R.vf_m[t+1, 2, :, :])
        vf_u_e_i = LinearInterpolation((P.grid_a, P.grid_b), R.vf_u[t+1, 1, :, :])
        vf_u_u_i = LinearInterpolation((P.grid_a, P.grid_b), R.vf_u[t+1, 2, :, :])

        # cycle over state variables
        for (i_y, y) = enumerate(P.grid_y), (i_a, a) = enumerate(P.grid_a), (i_b, b) = enumerate(P.grid_b)
        
            budget = y + a + P.δ * b * (1 + P.r)

            # unmatched agents
            function value_function_u(a_p, b_p) 
                if (a_p < P.min_a) | (a_p > P.max_a) | (b_p < P.min_b) | (b_p > P.max_b)
                    return -1/eps()
                else
                    # invert lt bond lom to get lt bond purchases
                    b_tilde = b_p - (1-P.δ) * b * (1+P.r)

                    if b_tilde < 0
                        return -1/eps()
                    end

                    # consumption is budget net of cash for tomorrow and lt bond purchases
                    consumption = budget - a_p - b_tilde

                    # instanteous value
                    result = u(consumption, P.σ)

                    # continuation value
                    result += P.β * R.γ * P.Π_y[i_y, 1] * vf_m_e_i[a_p, b_p]
                    result += P.β * R.γ * P.Π_y[i_y, 2] * vf_m_u_i[a_p, b_p]
                    result += P.β * (1 - R.γ) * P.Π_y[i_y, 1] * vf_u_e_i[a_p, b_p]
                    result += P.β * (1 - R.γ) * P.Π_y[i_y, 2] * vf_u_u_i[a_p, b_p]

                    return result
                end
            end

            # negative value function for minimizer
            function obj_u(x) 
                return -value_function_u(x[1], x[2])
            end

            # solve optimum
            optim = optimize(obj_u, [a, b], NelderMead())

            # save results
            R.vf_u[t, i_y, i_a, i_b] = -optim.minimum
            R.pf_u_a[t, i_y, i_a, i_b]  = optim.minimizer[1]
            R.pf_u_b_tilde[t, i_y, i_a, i_b]  = optim.minimizer[2] - (1-P.δ) * b * (1+P.r)
            R.pf_u_c[t, i_y, i_a, i_b] =  budget - R.pf_u_a[t, i_y, i_a, i_b] - R.pf_u_b_tilde[t, i_y, i_a, i_b]


            # figure out price
            function price_function(p, a_p, b_p, b_hat)
                if (p < 0)  | (b_p + b_hat > P.max_b)
                    return -1
                else
                    if (a_p + p + 0.001 >= P.max_a) & (a_p + p < P.max_a)
                        epsilon = -0.001
                    elseif (a_p + p >= P.max_a)
                        return -1
                    else 
                        epsilon = 0.001
                    end

                    # get value for intermediary
                    W_b_hat = W(b_hat, P.δ, R.β_I)

                    # compute value function at a' + p and b'
                    V_1 = R.γ * P.Π_y[i_y, 1] * vf_m_e_i[a_p + p, b_p]
                    V_1 += R.γ * P.Π_y[i_y, 2] * vf_m_u_i[a_p + p, b_p]
                    V_1 += (1 - R.γ) * P.Π_y[i_y, 1] * vf_u_e_i[a_p + p, b_p]
                    V_1 += (1 - R.γ) * P.Π_y[i_y, 2] * vf_u_u_i[a_p + p, b_p]

                    # compute value function at a' + p + epsilon and b'
                    V_1_eps = R.γ * P.Π_y[i_y, 1] * vf_m_e_i[a_p + p + epsilon, b_p]
                    V_1_eps += R.γ * P.Π_y[i_y, 2] * vf_m_u_i[a_p + p + epsilon, b_p]
                    V_1_eps += (1 - R.γ) * P.Π_y[i_y, 1] * vf_u_e_i[a_p + p + epsilon, b_p]
                    V_1_eps += (1 - R.γ) * P.Π_y[i_y, 2] * vf_u_u_i[a_p + p + epsilon, b_p]

                    # get numerical derivative at a' + p and b'
                    V_1_p = (V_1_eps - V_1)/epsilon

                    # compute value function at a' and b' + b_hat
                    V_2 = R.γ * P.Π_y[i_y, 1] * vf_m_e_i[a_p, b_p + b_hat]
                    V_2 += R.γ * P.Π_y[i_y, 2] * vf_m_u_i[a_p, b_p + b_hat]
                    V_2 += (1 - R.γ) * P.Π_y[i_y, 1] * vf_u_e_i[a_p, b_p + b_hat]
                    V_2 += (1 - R.γ) * P.Π_y[i_y, 2] * vf_u_u_i[a_p, b_p + b_hat]

                    return P.θ * V_1_p * (b_hat - p) - (1-P.θ) * (V_1 - V_2)                        
                end
            end


            # matched agents
            function value_function_m(a_p, l) 
                if (a_p < P.min_a) | (a_p > P.max_a) | (l < 0.0) | (l > 1.0)
                    return -1/eps()
                else

                    # get unsold, unmatured LT bonds
                    b_p = (1-l)*(1-P.δ)*b*(1+P.r)

                    # get liquidated amount of LT bonds
                    b_hat = l*(1-P.δ)*b*(1+P.r)

                    if b_hat > 0.0
                        p_temp = find_zeros(p -> price_function(p, a_p, b_p, b_hat), (P.min_b, max(1, b_hat)))
                        if length(p_temp) == 0
                            p = 0.0
                        else
                            p = p_temp[1]
                        end
                    else
                        p = 0.0
                    end

                    # consumption is budget net of cash for tomorrow
                    consumption = budget - a_p

                    # instanteous value
                    result = u(consumption, P.σ)

                    # continuation value
                    result += P.β * R.γ * P.Π_y[i_y, 1] * vf_m_e_i[a_p + p, b_p]
                    result += P.β * R.γ * P.Π_y[i_y, 2] * vf_m_u_i[a_p + p, b_p]
                    result += P.β * (1 - R.γ) * P.Π_y[i_y, 1] * vf_u_e_i[a_p + p, b_p]
                    result += P.β * (1 - R.γ) * P.Π_y[i_y, 2] * vf_u_u_i[a_p + p, b_p]

                    return result
                end
            end

            # negative value function for minimizer
            function obj_m(x) 
                return -value_function_m(x[1], x[2])
            end

            # solve optimum
            optim = optimize(obj_m, [a, 0.5], NelderMead())

            # if better to buy LT bonds
            if -optim.minimum < R.vf_u[t, i_y, i_a, i_b]
                R.vf_m[t, i_y, i_a, i_b] = R.vf_u[t, i_y, i_a, i_b]
                R.pf_m_a[t, i_y, i_a, i_b] = R.pf_u_a[t, i_y, i_a, i_b]
                R.pf_m_b_tilde[t, i_y, i_a, i_b] = R.pf_u_b_tilde[t, i_y, i_a, i_b]
                R.pf_m_c[t, i_y, i_a, i_b] =  R.pf_u_c[t, i_y, i_a, i_b]
                
            else
                a_p = optim.minimizer[1]
                l = optim.minimizer[2]
                b_hat = l*(1-P.δ)*b*(1+P.r)

                R.vf_m[t, i_y, i_a, i_b] = -optim.minimum
                R.pf_m_a[t, i_y, i_a, i_b] = a_p
                R.pf_m_l[t, i_y, i_a, i_b] = l
                R.pf_m_c[t, i_y, i_a, i_b] = budget - a_p
                R.pf_m_b_hat[t, i_y, i_a, i_b] =  b_hat

                p_temp = find_zeros(p -> price_function(p, a_p, (1-l)*(1-P.δ)*b*(1+P.r), b_hat), (P.min_b, max(1, l*(1-P.δ)*b*(1+P.r))))

                if length(p_temp) == 0
                    R.price[t, i_y, i_a, i_b] = 0.0
                else
                    R.price[t, i_y, i_a, i_b] = p_temp[1]
                end
            end
        end
    end
end


function Solve!(R::Results)
    Solve_terminal_period!(R)
    Solve_nonterminal_period!(R)
end