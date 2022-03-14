# Alex von Hafften
# ECON 810A - Project
# Bewley with long-term illiquid Assets

cd("/Users/vonhafften/Documents/UW Madison/problem_sets/econ_810a/project/")

include("model.jl")

using Parameters, Interpolations, LinearAlgebra, StatsBase

# Results
mutable struct Simulation_Results_2
    N::Int64                   # Number of trials
    y::Array{Float64}          # labor income
    m::Array{Int64}            # matched or unmatched.
    c::Array{Float64}          # consumption
    a::Array{Float64}          # cash holdings
    b::Array{Float64}          # LT bond holdings

end

# initialize simuation
function Initialize_Simulation(N::Int64)
    P = Primitives()

    y = zeros(N, P.T)
    m = fill(0, N, P.T)
    c = zeros(N, P.T)
    a = zeros(N, P.T)
    b = zeros(N, P.T)

    Simulation_Results_2(N, y, m, c, a, b)
end

function Simulate!(S::Simulation_Results_2, R::Results)
    P = Primitives()

    # first period. everyone starts with uniform distribution of cash and LT bonds on support
    S.a[:, 1] = rand(S.N) * (P.min_a + P.max_a)/2
    S.b[:, 1] = rand(S.N) * (P.min_b + P.max_b)/2

    # everyone starts employed
    S.y[:, 1] .= sample([1.0, 0.5], Weights(diag((P.Π_y)^1000)), S.N)

    # Interpolations
    pf_m_c_e_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_m_c[1, 1, :, :])
    pf_u_c_e_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_u_c[1, 1, :, :])
    pf_m_a_e_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_m_a[1, 1, :, :])
    pf_u_a_e_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_u_a[1, 1, :, :])
    price_e = LinearInterpolation((P.grid_a, P.grid_b), R.price[1, 1, :, :])
    pf_m_b_hat_e_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_m_b_hat[1, 1, :, :])
    pf_m_l_e_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_m_l[1, 1, :, :])
    pf_m_b_tilde_e_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_m_b_tilde[1, 1, :, :])
    pf_u_b_tilde_e_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_u_b_tilde[1, 1, :, :])

    # for first period
    for i = 1:S.N
        # if matched
        if rand() < P.γ
            S.m[i, 1] = 1
            S.c[i, 1] = pf_m_c_e_i[S.a[i, 1], S.b[i, 1]]

            # if choose to liquidate LT bonds
            if (pf_m_l_e_i[S.a[i, 1], S.b[i, 1]][1] > 0.0)
                S.a[i, 2] = pf_m_a_e_i[S.a[i, 1], S.b[i, 1]] + price_e[S.a[i, 1], S.b[i, 1]]
                S.b[i, 2] = (1- pf_m_l_e_i[S.a[i, 1], S.b[i, 1]]) * (1+P.r) * (1-P.δ) * S.b[i, 1]
            else # if to buy LT bonds
                S.a[i, 2] = pf_m_a_e_i[S.a[i, 1], S.b[i, 1]]
                S.b[i, 2] = (1+P.r) * (1-P.δ) * S.b[i, 1][1] + pf_m_b_tilde_e_i[S.a[i, 1], S.b[i, 1]]
            end
        else # if unmatched
            S.m[i, 1] = 0
            S.c[i, 1] = pf_u_c_e_i[S.a[i, 1], S.b[i, 1]]
            S.a[i, 2] = pf_u_a_e_i[S.a[i, 1], S.b[i, 1]]
            S.b[i, 2] = (1+P.r) * (1-P.δ) * S.b[i, 1][1] + pf_u_b_tilde_e_i[S.a[i, 1], S.b[i, 1]]
        end
    end

    # other periods
    for t = 2:P.T
        pf_m_c_e_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_m_c[t, 1, :, :])
        pf_m_c_u_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_m_c[t, 2, :, :])

        pf_u_c_e_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_u_c[t, 1, :, :])
        pf_u_c_u_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_u_c[t, 2, :, :])

        pf_m_a_e_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_m_a[t, 1, :, :])
        pf_m_a_u_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_m_a[t, 2, :, :])

        pf_u_a_e_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_u_a[t, 1, :, :])
        pf_u_a_u_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_u_a[t, 2, :, :])

        pf_m_b_tilde_e_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_m_b_tilde[t, 1, :, :])
        pf_m_b_tilde_u_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_m_b_tilde[t, 2, :, :])

        pf_u_b_tilde_e_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_u_b_tilde[t, 1, :, :])
        pf_u_b_tilde_u_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_u_b_tilde[t, 2, :, :])

        pf_m_b_hat_e_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_m_b_hat[t, 1, :, :])
        pf_m_b_hat_u_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_m_b_hat[t, 2, :, :])

        pf_m_l_e_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_m_b_hat[t, 1, :, :])
        pf_m_l_u_i = LinearInterpolation((P.grid_a, P.grid_b), R.pf_m_b_hat[t, 2, :, :])

        price_e = LinearInterpolation((P.grid_a, P.grid_b), R.price[t, 1, :, :])
        price_u = LinearInterpolation((P.grid_a, P.grid_b), R.price[t, 2, :, :])

        for i = 1:S.N
            # draw employment shock
            if S.y[i, t-1] == 1.0
                if rand() < P.Π_y[1, 1]
                    S.y[i, t] = 1.0
                else
                    S.y[i, t] = 0.5
                end
            else
                if rand() < P.Π_y[2, 1]
                    S.y[i, t] = 1.0
                else
                    S.y[i, t] = 0.5
                end
            end

            # if employed
            if S.y[i, t] == 1.0
                # if matched
                if rand() < P.γ
                    S.m[i, t] = 1
                    S.c[i, t] = pf_m_c_e_i[S.a[i, t], S.b[i, t]]
                    
                    if t < P.T
                        # if choose to liquidate LT bonds
                        if (pf_m_l_e_i[S.a[i, t], S.b[i, t]] > 0.0)
                            S.a[i, t+1] = max(pf_m_a_e_i[S.a[i, t], S.b[i, t]] + price_e[S.a[i, t], S.b[i, t]], 0.0)
                            S.b[i, t+1] = max((1- pf_m_l_e_i[S.a[i, t], S.b[i, t]][1]) * (1+P.r) * (1-P.δ) * S.b[i, t], 0.0)

                        else # if to buy LT bonds
                            S.a[i, t+1] = max(pf_m_a_e_i[S.a[i, t], S.b[i, t]][1], 0.0)
                            S.b[i, t+1] = max( (1+P.r) * (1-P.δ) * S.b[i, t] + pf_m_b_tilde_e_i[S.a[i, t], S.b[i, t]], 0.0)
                        end
                    end

                else # if unmatched
                    S.m[i, t] = 0
                    S.c[i, t] = pf_u_c_e_i[S.a[i,t], S.b[i,t]]
                    
                    if t < P.T
                        S.a[i, t + 1] = max(pf_u_a_e_i[S.a[i, t], S.b[i, t]], 0.0)
                        S.b[i, t + 1] = max((1+P.r) * (1-P.δ) * S.b[i, t] + pf_u_b_tilde_e_i[S.a[i, t], S.b[i, t]], 0.0)
                    end
                end
            end

            # if unemployed
            if S.y[i, t] == 0.5
                # if matched
                if rand() < P.γ
                    S.m[i, t] = 1
                    S.c[i, t] = pf_m_c_u_i[S.a[i, t], S.b[i, t]]

                    if t < P.T
        
                        # if choose to liquidate LT bonds
                        if (pf_m_l_u_i[S.a[i, t], S.b[i, t]] > 0.0)
                            S.a[i, t+1] = max(pf_m_a_u_i[S.a[i, t], S.b[i, t]][1] + price_u[S.a[i, t], S.b[i, t]], 0.0)
                            S.b[i, t+1] = max((1- pf_m_l_u_i[S.a[i, t], S.b[i, t]][1]) * (1+P.r) * (1-P.δ) * S.b[i, t], 0.0)

                        else # if to buy LT bonds
                            S.a[i, t+1] = max(pf_m_a_u_i[S.a[i, t], S.b[i, t]], 0.0)
                            S.b[i, t+1] = max((1+P.r) * (1-P.δ) * S.b[i, t] + pf_m_b_tilde_u_i[S.a[i, t], S.b[i, t]], 0.0)
                        end
                    end

                else # if unmatched
                    S.m[i, t] = 0
                    S.c[i, t] = pf_u_c_u_i[S.a[i, t], S.b[i, t]]

                    if t < P.T
                        S.a[i, t + 1] = max(pf_u_a_u_i[S.a[i, t], S.b[i, t]], 0.0)
                        S.b[i, t + 1] = max((1+P.r) * (1-P.δ) * S.b[i, t] + pf_u_b_tilde_u_i[S.a[i, t], S.b[i, t]], 0.0)
                    end
                end
            end
        end
    end
end