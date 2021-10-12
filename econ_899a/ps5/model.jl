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

# simulate business cycles
function simulate_Z()
    @unpack T, z, Π_z = Primitives()
    @assert size(Π_z)[1] == size(Π_z)[2] # square required

    # setup the simulation
    X = fill(0, T) # allocate memory
    X[1] = 0 # set the initial state
   
    # iterate through each period
    # X is recession flags
    # X[t] == 0 in boom
    # X[t] == 1 in recession
    # rand() is uniform [0, 1), so if it is greater than Π_z[X[i-1], 1]
    # then this period is second state.
    for t in 2:T  
        X[t] = Π_z[X[t-1]+1, 1] < rand()
    end
    
    return z[X .+ 1]
end

# simulate employment based on business cycles
function simulate_E(Z)
    @unpack T, N, ε, Π_ε, Π_ε_star = Primitives()
    @assert size(Π_ε)[1] == size(Π_ε)[2]

    # X_Z is recession flags
    # X_Z[t] == 0 in boom
    # X_Z[t] == 1 in recession
    X_Z = Int64.(Z .== z[2])

    # setup the simulation
    X = fill(0, T, N) # allocate memory

    # iterate through each agent
    # X is employed flag
    # X[t,i] == 0 in unemployed
    # X[t,i] == 1 in employed
    for i in 1:N
        
        # Choose initial employment status based on stationary distribution in initial period.
        X[1, i] = Int64(Π_ε_star[X_Z[1]+ 1]/(Π_ε_star[X_Z[1]+ 1] + Π_ε_star[X_Z[1]+3]) > rand())

        # iterate through future periods
        for t in 2:T
            
            X_Z[t]
            X[t, i]

        end
    end
    return X
end

function Initialize()
    Z   = simulate_Z()
    E   = simulate_E(Z)
    a_0 = 0.095
    a_1 = 0.999
    b_0 = 0.085
    b_1 = 0.999

end




