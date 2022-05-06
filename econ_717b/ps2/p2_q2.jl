# Alex von Hafften
# ECON 717
# May 2, 2022
# Part 2 - Question 2

using Parameters, Distributions, StatsBase

cd("/Users/vonhafften/Documents/UW Madison/problem_sets/econ_717b/ps2/")

@with_kw struct paramters
    σ_00 = 1.0
    σ_11 = 1.0
    σ_vv = 1.0
    σ_01 = 0.5
    σ_0v = 0.3
    σ_1v = 0.7
    Σ    = [σ_00 σ_01 σ_0v; σ_01 σ_11 σ_1v; σ_0v σ_1v σ_vv]
    π    = 0.3 # pr(z = 1)
end

# data object
mutable struct data
    # number of trials
    N::Int64              
    
    # unobservables
    u_0::Vector{Float64}
    u_1::Vector{Float64}
    u_v::Vector{Float64}
    y_0::Vector{Float64}
    y_1::Vector{Float64}
    v::Vector{Float64}
    is_complier::BitVector

    # observables
    d::Vector{Int64}
    z::Vector{Float64}
    y::Vector{Float64}
end

# simulates data based on the seed
function simulate_data(; seed = -1, N = 10000, part_d = false)
    p = paramters()

    # draw unobservables
    u_distribution = MvNormal([0, 0, 0], p.Σ)
    u_draws = rand(u_distribution, N)'
    u_0 = u_draws[:, 1]
    u_1 = u_draws[:, 2]
    u_v = u_draws[:, 3]

    # draw instruments
    z = sample([0, 1], Weights([1.0 - p.π, p.π]), N)
    
    # potential outcomes
    y_0 = 1.0 .+ u_0
    y_1 = 4.0 .+ u_1

    if part_d
        v = -1.0 .+ 2.0 .* z .+ u_v
        v_0 = -1.0 .+ u_v
        v_1 = -1.0 .+ 2.0 .+ u_v
    else
        v = -1.0 .+ 3.0 .* z .+ u_v
        v_0 =  -1.0 .+ u_v
        v_1 =  -1.0 .+ 3.0 .+ u_v
    end

    is_complier = (v_0 .< 0.0) .* (v_1 .> 0.0)
    d = v .>= 0
    y = y_1 .* d + y_0 .* (1 .- d)

    data(N, u_0, u_1, u_v, y_0, y_1, v, is_complier, d, z, y)
end

##################################################
# part a - generate data
##################################################

d = simulate_data()

##################################################
# part b - ate, atet, ateu, ols/naive, 
#          direct/reduced reform/itt, iv
##################################################

# ate
mean(d.y_1 .- d.y_0)

# atet
mean(d.y_1[d.v .> 0] .- d.y_0[d.v .> 0])

# ateu
mean(d.y_1[d.v .< 0] .- d.y_0[d.v .< 0])

# ols/naive estimator
mean(d.y[d.v .> 0]) - mean(d.y[d.v .< 0])

# direct/reduced form/itt estimator
mean(d.y[d.z .== 1]) - mean(d.y[d.z .== 0])

# iv estimator
(mean(d.y[d.z .== 1]) - mean(d.y[d.z .== 0])) / (mean(d.d[d.z .== 1]) - mean(d.d[d.z .== 0]))

##################################################
# part c - what fraction are compliers?
##################################################

sum(d.is_complier) ./ d.N

##################################################
# part d - different v function
##################################################

d = simulate_data(; part_d = true)

# ate
mean(d.y_1 .- d.y_0)

# atet
mean(d.y_1[d.v .> 0] .- d.y_0[d.v .> 0])

# ateu
mean(d.y_1[d.v .< 0] .- d.y_0[d.v .< 0])

# ols/naive estimator
mean(d.y[d.v .> 0]) - mean(d.y[d.v .< 0])

# direct/reduced form/itt estimator
mean(d.y[d.z .== 1]) - mean(d.y[d.z .== 0])

# iv estimator
(mean(d.y[d.z .== 1]) - mean(d.y[d.z .== 0])) / (mean(d.d[d.z .== 1]) - mean(d.d[d.z .== 0]))

# fraction of compliers
sum(d.is_complier) ./ d.N

