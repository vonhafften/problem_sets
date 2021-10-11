# Krusell and Smith (1998)
# Alex von Hafften
# October 13, 2021

# ECON 899A Computational Economics
# Problem set 5

# This file contains the functions to estimate the model.
# ./run.jl calls this file and runs the model.

using Parameters

include("transition_functions.jl");

# Parameters
@with_kw struct Primitives
    β::Float64               = 0.99                    # discount rate
    α::Float64               = 0.36                    # elasticity of capital in production
    z::Array{Float64, 1}     = [1.01, 0.99]            # TFP states
    Π_z::Array{Float64}      = get_Π_z()               # aggregate shock transition matrix
    δ::Float64               = 0.025                   # Capital depreciation
    ε::Array{Float64, 1}     = [0, 1]                  # Idiosyncratic employment opportunities
    e::Float64               = 0.3271                  # Labor efficiency per unit of time worked
    Π_ε::Array{Float64}      = get_Π_ε()               # employment transition matrix
    Π_ε_star::Array{Float64} = compute_Π_star(Π_ε)     # stationary distribution of Π
    T::Int64                 = 11000                   # Number of periods in simulation
    N::Int64                 = 5000                    # Number of individuals in simulation
end

mutable struct Results
    Z::Array{Float64, 1}
    E::Array{Float64, 2}
    a_0::Float64
    a_1::Float64
    b_0::Float64
    b_1::Float64
end

# simulate business cycle
function simulate_Z(z, Π_z, T)

    @assert size(Π_z)[1] == size(Π_z)[2] # square required
    N = size(Π_z)[1] # should be square

    # setup the simulation
    X = fill(0, T) # allocate memory, or zeros(Int64, sample_size)
    X[1] = 1 # set the initial state

    # iterate through each period
    for i in 2:T
        
        # rand() is uniform [0, 1), so if less than Π_z[X[i-1], 1], then this period is first state.
        if Π_z[X[i-1], 1] > rand()
            X[i] = 1
        else
            X[i] = 2
        end
    end
    
    return z[X]
end



function Initialize()
    @unpack T, N, z, ε, Π_z, Π_ε = Primitives()

    Z   = simulate_Z(z, Π_z, T)
    E   = simulate_E(ε, Π_z_ε, T, N)
    a_0 = 

end




