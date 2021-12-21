# Problem Set 3 - BLP
# Computational Economics
# Taught by JF Houde

# Alex von Hafften

using DataFrames, StatFiles, LinearAlgebra, Statistics

# function to get data for particular market
function get_segment_data(market)
    df = filter(row -> row.Year == market, df_characteristics)

    δ_iia = df.delta_iia
    S     = df.share
    P     = df.price
    Y     = df_types

    return δ_iia, S, Y, P
end

# function to invert demand
# tolerence[1] is when to stop
# tolerence[2] is when to switch from contraction mapping to Newton's method
function invert_demand(market, λ_p::Float64; tolerence = [1e-12, 1])
    
    δ_iia, S, Y, P = get_segment_data(market)
    R = length(Y)
    J = length(S)

    # compute idiosyncratic utility
    μ = λ_p * P * Y'

    err, err_list = 100, []
    i, maxiter = 1, 1000

    δ_0 = copy(δ_iia)
    δ_1 = copy(δ_iia)

    while (err > tolerence[1]) & (i < maxiter)
    
        # computes choice probability
        Λ = exp.(δ_0 .+ μ)
        Σ = Λ ./ sum(1 .+ Λ; dims = 1)
        σ = sum(Σ; dims =2) / R

        # updates guess
        # if err is large use contraction mapping
        if (err > tolerence[2])
            δ_1 = δ_0 + log.(S) - log.(σ)
        else # if err is small use Newton
            # compute jacobian
            Δ = (1 / R) * ((I(J) .* (Σ * (1 .- Σ)')) - ((1 .- I(J)) .* (Σ * Σ'))) ./ σ
            δ_1 = δ_0 + inv(Δ) * ( log.(S) - log.(σ))
        end

        err = norm(δ_0 - δ_1)
        push!(err_list, err)

        δ_0 = copy(δ_1)

        i += 1
    
    end
    return δ_0, err_list
end

function invert_demand(λ_p::Float64)
    δ = []
    
    for m in markets
        δ = vcat(δ, invert_demand(m, λ_p)[1])
    end

    return δ
end

function compute_ρ(λ_p::Float64, W::Array{Float64})
    δ = invert_demand(λ_p)

    β = inv((X'*Z) *W*(Z'X)) *(X'*Z) *W*Z'*δ

    return δ - X*β
end

# compute gmm objective function
function gmm_obj(λ_p::Float64, W::Array{Float64})
    ρ = compute_ρ(λ_p, W)

    return (ρ'*Z*W*Z'*ρ)[1,1]

end

