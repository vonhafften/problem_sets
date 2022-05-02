# Alex von Hafften
# ECON 717
# May 2, 2022

using Parameters, Distributions

cd("/Users/vonhafften/Documents/UW Madison/problem_sets/econ_717b/ps2/")

# data object
mutable struct data
    N::Int64              # number of trials

    # unobservables
    a::Vector{Float64}    # ability
    ε::Vector{Float64}    # shock to earnings
    η::Vector{Float64}    # shock to selection
    
    # observables
    ln_w::Vector{Float64} # wages
    s::Vector{Float64}    # schooling
    z_1::Vector{Float64}  # cost of schooling 1
    z_2::Vector{Float64}  # cost of schooling 2
    z_3::Vector{Float64}  # cost of schooling 3
end

# simulates data based on the seed
function simulate_data(; seed = -1, N = 2000)

    # draw random variables
    ε   =  rand(Normal(0.0,  0.5), N)
    a   =  rand(Normal(0.0,  4.0), N)
    η   =  rand(Normal(0.0,  1.0), N)
    z_1 =  rand(Normal(0.0,  0.1), N)
    z_2 =  rand(Normal(0.0, 25.0), N)
    z_3 = rand(Uniform(0.0,  1.0), N)

    # schooling and ln w
    s    = 3.0 .* a .+ z_1 .+ z_2 .+ η
    ln_w = 1 .+ 0.05 .* s + 0.1 .* a .+ ε 

    data(N, a, ε, η, ln_w, s, z_1, z_2, z_3)
end


# part a - generate data
data_1 = simulate_data()

# part b - ols
Y = data_1.ln_w
X = hcat(ones(data_1.N), data_1.s)

β_OLS = inv(X' * X) * X' * Y

# part c - see write-up

# part d - 2sls
Z = hcat(ones(data_1.N), data_1.z_1, data_1.z_2, data_1.z_2)
