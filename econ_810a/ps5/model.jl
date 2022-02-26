# Alex von Hafften
# ECON 810: Advanced Macro
# PS 5 - Code for solving model
# Professor Carter Braxton

using Parameters, DataFrames, QuantEcon, CSV, LinearAlgebra

cd("/Users/vonhafften/Documents/uw_madison/problem_sets/econ_810a/ps5/")

# Primitive structure
@with_kw struct Primitives

    # parameters
    R_f::Float64 = 1.04^6   # gross risk-free rate
    β::Float64   = 1/R_f    # discount rate
    ρ_h::Float64 = (0.97)^6 # persistence of human capital shocks
    σ_h::Float64 = 0.13     # variance of human capital shocks
    ω_c::Float64 = 0.5      # relative weight of previous hc and current investment in child's human capital evolution
    γ::Float64   = 1/10     # dollars to human capital exchange rate
    θ::Float64   = 0.5      # altruism
    σ::Float64   = 2.0      # coefficient of relative risk aversion

    # Age profile of earnings
    κ::Vector{Float64} = CSV.read("age_profile.csv", DataFrame).kappa

    # human capital grid
    N_h::Int64              = 11                      # number of points
    tauchen_h               = tauchen(N_h, ρ_h, σ_h) # tauchen object
    Π_h::Matrix{Float64}    = tauchen_h.p            # transition probabilities
    grid_h::Vector{Float64} = tauchen_h.state_values # grid for human capital
    min_h::Float64          = minimum(grid_h)        # lower bound
    max_h::Float64          = maximum(grid_h)        # upper bound
    Π_h_0                   = vcat(ones(fld(N_h, 2)), zeros(cld(N_h, 2)))/fld(N_h, 2) # child hc is drawn from uniform distribution on bottom half of distribution

    # bond grid
    min_b::Float64            = 0.0                               # lower bound
    max_b::Float64            = 100000.0                          # upper bound
    N_b::Int64                = 101                               # number of points
    grid_srl_b                = range(min_b, max_b; length = N_b) # step range
    grid_b::Array{Float64, 1} = collect(grid_srl_b)               # grid array

end

# results
mutable struct Results

    # value functions
    vf_4::Array{Float64}
    vf_5::Array{Float64}
    vf_6::Array{Float64}
    vf_7::Array{Float64}
    vf_8::Array{Float64}
    vf_9::Array{Float64}
    vf_10::Array{Float64}
    vf_11::Array{Float64}
    vf_12::Array{Float64}

    # policy functions
    # assets
    pf_b_4::Array{Float64}
    pf_b_5::Array{Float64}
    pf_b_6::Array{Float64}
    pf_b_7::Array{Float64}
    pf_b_8::Array{Float64}
    pf_b_9::Array{Float64}
    pf_b_10::Array{Float64}
    pf_b_11::Array{Float64}

    # investment in child
    pf_i_5::Array{Float64}
    pf_i_6::Array{Float64}
    pf_i_7::Array{Float64}
    pf_i_8::Array{Float64}

    # transfer to child
    pf_τ_9::Array{Float64}
end

function Initialize()
    P = Primitives()

    # value functions
    vf_4 = zeros(P.N_b, P.N_h)
    vf_5 = zeros(P.N_b, P.N_h, P.N_h)
    vf_6 = zeros(P.N_b, P.N_h, P.N_h)
    vf_7 = zeros(P.N_b, P.N_h, P.N_h)
    vf_8 = zeros(P.N_b, P.N_h, P.N_h)
    vf_9 = zeros(P.N_b, P.N_h, P.N_h)
    vf_10 = zeros(P.N_b, P.N_h)
    vf_11 = zeros(P.N_b, P.N_h)
    vf_12 = zeros(P.N_b, P.N_h)

    # policy functions
    # assets
    pf_b_4 = zeros(P.N_b, P.N_h)
    pf_b_5 = zeros(P.N_b, P.N_h, P.N_h)
    pf_b_6 = zeros(P.N_b, P.N_h, P.N_h)
    pf_b_7 = zeros(P.N_b, P.N_h, P.N_h)
    pf_b_8 = zeros(P.N_b, P.N_h, P.N_h)
    pf_b_9 = zeros(P.N_b, P.N_h, P.N_h)
    pf_b_10 = zeros(P.N_b, P.N_h)
    pf_b_11 = zeros(P.N_b, P.N_h)

    # investment in child
    pf_i_5 = zeros(P.N_b, P.N_h, P.N_h)
    pf_i_6 = zeros(P.N_b, P.N_h, P.N_h)
    pf_i_7 = zeros(P.N_b, P.N_h, P.N_h)
    pf_i_8 = zeros(P.N_b, P.N_h, P.N_h)

    # transfer to child
    pf_τ_9 = zeros(P.N_b, P.N_h, P.N_h)

    Results(vf_4, vf_5, vf_6, vf_7, vf_8, vf_9, vf_10, vf_11, vf_12, 
            pf_b_4, pf_b_5, pf_b_6, pf_b_7, pf_b_8, pf_b_9, pf_b_10, pf_b_11,
            pf_i_5, pf_i_6, pf_i_7, pf_i_8,
            pf_τ_9)
end

# solve problems

# utility function
function u(c::Float64, σ::Float64)
    if (c > 0)
        return (c^(1-σ))/(1-σ)
    else
        return -1/eps()
    end
end

# solve t = 12 problem
function solve_12!(R::Results)
    P = Primitives()

    for (i_b, b) = enumerate(P.grid_b), (i_h, h) = enumerate(P.grid_h)
        budget = exp(P.κ[12-3] + h) + b*P.R_f
        R.vf_12[i_b, i_h] = u(budget, P.σ)
    end
    
end

# solve t = 11 problem
function solve_11!(R::Results)
    P = Primitives()

    for (i_b, b) = enumerate(P.grid_b), (i_h, h) = enumerate(P.grid_h)
        budget = exp(P.κ[11-3] + h) + b*P.R_f

        candidate_max = -1/eps()
        
        for (i_b_p, b_p) = enumerate(P.grid_b)
            consumption = budget - b_p

            # instanteous value
            value = u(consumption, P.σ)

            # continuation value
            value += P.β * R.vf_12[i_b_p, :]' * P.Π_h[i_h, :]

            # save if new candidate
            if value > candidate_max
                R.vf_11[i_b, i_h] = value
                R.pf_b_11[i_b, i_h] = b_p
                candidate_max = value
            end
        end
    end
end

# solve t = 10 problem
function solve_10!(R::Results)
    P = Primitives()

    for (i_b, b) = enumerate(P.grid_b), (i_h, h) = enumerate(P.grid_h)
        budget = exp(P.κ[10-3] + h) + b*P.R_f

        candidate_max = -1/eps()
        
        for (i_b_p, b_p) = enumerate(P.grid_b)
            consumption = budget - b_p

            # instanteous value
            value = u(consumption, P.σ)

            # continuation value
            value += P.β * R.vf_11[i_b_p, :]' * P.Π_h[i_h, :]

            # save if new candidate
            if value > candidate_max
                R.vf_10[i_b, i_h] = value
                R.pf_b_10[i_b, i_h] = b_p
                candidate_max = value
            end
        end
    end
end

# solve t = 9 problem
function solve_9!(R::Results)
    P = Primitives()

    # cycle over assets today, human capital today, and child's human capital today
    for (i_b, b) = enumerate(P.grid_b), (i_h, h) = enumerate(P.grid_h), (i_h_c, h_c) = enumerate(P.grid_h)
        budget = exp(P.κ[9-3] + h) + b*P.R_f

        candidate_max = -1/eps()
        
        # cycle over assets tomorrow and transfers to child
        for (i_b_p, b_p) = enumerate(P.grid_b), (i_τ, τ) = enumerate(P.grid_b)

            # trying to speed things up
            if b_p + τ > budget
                continue
            end

            consumption = budget - b_p - τ

            # instanteous value
            value = u(consumption, P.σ)

            # utility from transfer to child
            value += P.θ * R.vf_4[i_τ, i_h_c]

            # continuation value
            value += P.β * R.vf_10[i_b_p, :]' * P.Π_h[i_h, :]

            # save if new candidate
            if value > candidate_max
                R.vf_9[i_b, i_h, i_h_c] = value
                R.pf_b_9[i_b, i_h, i_h_c] = b_p
                R.pf_τ_9[i_b, i_h, i_h_c] = τ
                candidate_max = value
            end
        end
    end
end


# solve t = 8 problem
function solve_8!(R::Results)
    P = Primitives()

    # cycle over assets today, human capital today, and child's human capital today
    for (i_b, b) = enumerate(P.grid_b), (i_h, h) = enumerate(P.grid_h), (i_h_c, h_c) = enumerate(P.grid_h)
        budget = exp(P.κ[8-3] + h) + b*P.R_f

        candidate_max = -1/eps()
        
        # cycle over assets tomorrow and child's human capital
        for (i_b_p, b_p) = enumerate(P.grid_b), (i_h_c_p, h_c_p) =  enumerate(P.grid_h)

            # investment required to get child human capital to i_h_c_p
            investment = (exp(h_c_p) - (1 - P.ω_c) * exp(h_c)) / (P.γ * P.ω_c)

            # prevents negative investment in childs human capital
            if investment < 0
                continue
            end

            consumption = budget - b_p - investment

            # instanteous value
            value = u(consumption, P.σ)

            # continuation value
            value += P.β * R.vf_9[i_b_p, :, i_h_c_p]' * P.Π_h[i_h, :]

            # save if new candidate
            if value > candidate_max
                R.vf_8[i_b, i_h, i_h_c] = value
                R.pf_b_8[i_b, i_h, i_h_c] = b_p
                R.pf_i_8[i_b, i_h, i_h_c] = investment
                candidate_max = value
            end
        end
    end
end


# solve t = 7 problem
function solve_7!(R::Results)
    P = Primitives()

    # cycle over assets today, human capital today, and child's human capital today
    for (i_b, b) = enumerate(P.grid_b), (i_h, h) = enumerate(P.grid_h), (i_h_c, h_c) = enumerate(P.grid_h)
        budget = exp(P.κ[7-3] + h) + b*P.R_f

        candidate_max = -1/eps()

        # cycle over assets tomorrow and transfers to child
        for (i_b_p, b_p) = enumerate(P.grid_b), (i_h_c_p, h_c_p) = enumerate(P.grid_h)

            # investment required to get child human capital to i_h_c_p
            investment = (exp(h_c_p) - (1 - P.ω_c) * exp(h_c)) / (P.γ * P.ω_c)

            # prevents negative investment in childs human capital
            if investment < 0
                continue
            end

            consumption = budget - b_p - investment

            # instanteous value
            value = u(consumption, P.σ)

            # continuation value
            value += P.β * R.vf_8[i_b_p, :, i_h_c_p]' * P.Π_h[i_h, :]

            # save if new candidate
            if value > candidate_max
                R.vf_7[i_b, i_h, i_h_c] = value
                R.pf_b_7[i_b, i_h, i_h_c] = b_p
                R.pf_i_7[i_b, i_h, i_h_c] = investment
                candidate_max = value
            end
        end
    end
end


# solve t = 6 problem
function solve_6!(R::Results)
    P = Primitives()

    # cycle over assets today, human capital today, and child's human capital today
    for (i_b, b) = enumerate(P.grid_b), (i_h, h) = enumerate(P.grid_h), (i_h_c, h_c) = enumerate(P.grid_h)
        budget = exp(P.κ[6-3] + h) + b*P.R_f

        candidate_max = -1/eps()
        
        # cycle over assets tomorrow and transfers to child
        for (i_b_p, b_p) = enumerate(P.grid_b), (i_h_c_p, h_c_p) = enumerate(P.grid_h)

            # investment required to get child human capital to i_h_c_p
            investment = (exp(h_c_p) - (1 - P.ω_c) * exp(h_c)) / (P.γ * P.ω_c)

            # prevents negative investment in childs human capital
            if investment < 0
                continue
            end

            consumption = budget - b_p - investment

            # instanteous value
            value = u(consumption, P.σ)

            # continuation value
            value += P.β * R.vf_7[i_b_p, :, i_h_c_p]' * P.Π_h[i_h, :]

            # save if new candidate
            if value > candidate_max
                R.vf_6[i_b, i_h, i_h_c] = value
                R.pf_b_6[i_b, i_h, i_h_c] = b_p
                R.pf_i_6[i_b, i_h, i_h_c] = investment
                candidate_max = value
            end
        end
    end
end


# solve t = 5 problem
function solve_5!(R::Results)
    P = Primitives()

    # cycle over assets today, human capital today, and child's human capital today
    for (i_b, b) = enumerate(P.grid_b), (i_h, h) = enumerate(P.grid_h), (i_h_c, h_c) = enumerate(P.grid_h)
        budget = exp(P.κ[5-3] + h) + b*P.R_f

        candidate_max = -1/eps()
        
        # cycle over assets tomorrow and transfers to child
        for (i_b_p, b_p) = enumerate(P.grid_b), (i_h_c_p, h_c_p) = enumerate(P.grid_h)

            # investment required to get child human capital to i_h_c_p
            investment = (exp(h_c_p) - (1 - P.ω_c) * exp(h_c)) / (P.γ * P.ω_c)

            # prevents negative investment in childs human capital
            if investment < 0
                continue
            end

            consumption = budget - b_p - investment

            # instanteous value
            value = u(consumption, P.σ)

            # continuation value
            value += P.β * R.vf_6[i_b_p, :, i_h_c_p]' * P.Π_h[i_h, :]

            # prevents negative investment in childs human capital
            if investment < 0
                value = -1/eps()
            end

            # save if new candidate
            if value > candidate_max
                R.vf_5[i_b, i_h, i_h_c] = value
                R.pf_b_5[i_b, i_h, i_h_c] = b_p
                R.pf_i_5[i_b, i_h, i_h_c] = investment
                candidate_max = value
            end
        end
    end
end

# solve t = 4 problem
function solve_4!(R::Results)
    P = Primitives()

    # cycle over assets today, human capital today, and child's human capital today
    for (i_b, b) = enumerate(P.grid_b), (i_h, h) = enumerate(P.grid_h)
        budget = exp(P.κ[4-3] + h) + b*P.R_f

        candidate_max = -1/eps()
        
        # cycle over assets tomorrow and transfers to child
        for (i_b_p, b_p) = enumerate(P.grid_b)

            consumption = budget - b_p

            # instanteous value
            value = u(consumption, P.σ)

            # continuation value
            value += P.β * (R.vf_5[i_b_p, :, :] * P.Π_h[i_h, :])' * P.Π_h_0

            # save if new candidate
            if value > candidate_max
                R.vf_4[i_b, i_h] = value
                R.pf_b_4[i_b, i_h] = b_p
                candidate_max = value
            end
        end
    end
end

function Solve!(R::Results) 
    i = 1
    error = 100
    max_iter = 10
    tolerence = 1e-5

    while error > tolerence
        if i > max_iter
            break
        end
        
        println("Iteration #", i)

        # save previous vf_4
        vf_4_previous = copy(R.vf_4)

        # backward induction
        print("Working on t = 12... ")
        solve_12!(R)
        println("Done.")

        print("Working on t = 11... ")
        solve_11!(R)
        println("Done.")

        print("Working on t = 10... ")
        solve_10!(R)
        println("Done.")

        print("Working on t = 9... ")
        solve_9!(R)
        println("Done.")

        print("Working on t = 8... ")
        solve_8!(R)
        println("Done.")

        print("Working on t = 7... ")
        solve_7!(R)
        println("Done.")

        print("Working on t = 6... ")
        solve_6!(R)
        println("Done.")

        print("Working on t = 5... ")
        solve_5!(R)
        println("Done.")

        print("Working on t = 4... ")
        solve_4!(R)
        println("Done.")

        # evaluate sup norm
        error = maximum(abs.(vf_4_previous .- R.vf_4))

        println("Error: ", error)

        i += 1
    end
end