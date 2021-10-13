# Krusell and Smith (1998)
# Alex von Hafften
# October 13, 2021

# ECON 899A Computational Economics
# Problem set 5

# This file contains the functions to estimate the model.
# ./run.jl calls this file and runs the model.

using Parameters

include("transition_functions.jl");

####################################################################################################################
##################################### Parameter and result structures ##############################################
####################################################################################################################

@with_kw struct Primitives
    
    # Preference parameters
    β::Float64                = 0.99                    # discount rate

    # Production technology parameters
    α::Float64                = 0.36                    # elasticity of capital in production
    z::Array{Float64, 1}      = [1.01, 0.99]            # TFP states
    Π_z::Array{Float64}       = get_Π_z()               # aggregate shock transition matrix
    δ::Float64                = 0.025                   # Capital depreciation
    ε::Array{Float64, 1}      = [0, 1]                  # Idiosyncratic employment opportunities
    e::Float64                = 0.3271                  # Labor efficiency per unit of time worked
    Π_ε::Array{Float64}       = get_Π_ε()               # employment transition matrix
    Π_ε_star::Array{Float64}  = compute_Π_star(Π_ε)     # stationary distribution of Π

    # Simulation parameters
    T::Int64                  = 11000                   # Number of periods in simulation
    N::Int64                  = 5000                    # Number of individuals in simulation

    # Aggregate capital grid
    K_min::Float64                = 11.0                                   # Minimum of aggregate capital grid
    K_max::Float64                = 15.0                                   # Maximum of aggregate capital grid
    K_length::Int64               = 17                                     # Number of aggregate capital grid points
    K_grid_srl::Array{Float64, 1} = range(K_min, K_max; length = K_length) # Aggregate capital grid step range length
    K_grid::Array{Float64, 1}     = collect(K_grid_srl)                    # Aggregate capital grid

    # Individual capital grid
    k_min::Float64                = 0.0                                    # Minimum of individual capital grid
    k_max::Float64                = 15.0                                   # Maximum of individual capital grid
    k_length::Int64               = 61                                     # Number of individual capital grid points
    k_grid_srl::Array{Float64, 1} = range(K_min, K_max; length = K_length) # individual capital grid step range length
    k_grid::Array{Float64, 1}     = collect(K_grid_srl)                    # individual capital grid
end

mutable struct Results
    Z::Array{Float64, 1} # Aggregate productivity
    E::Array{Float64, 2} # Employment 
    a_0::Float64         # Intercept in boom capital regression
    a_1::Float64         # Slope in boom capital regression
    b_0::Float64         # Intercept in recession capital regression
    b_1::Float64         # Slope in recession capital regression
    K::Array{Float64, 1} # Aggregate capital
end

####################################################################################################################
####### Functions to simulate business cycles and employment based on business cycles. #############################
####################################################################################################################

# simulate business cycles
function simulate_Z()
    @unpack T, z, Π_z = Primitives()
    @assert size(Π_z)[1] == size(Π_z)[2] # square required

    # setup the simulation
    X = fill(0, T) # allocate memory
    X[1] = 0 # set the initial state
   
    # iterate through each period
    # X is recession indicator (X[t] == 0 in boom and X[t] == 1 in recession)
    for t in 2:T  
        X[t] = Π_z[X[t-1]+1, 1] < rand()
    end
    
    return z[X .+ 1]
end

# simulate employment based on business cycles
function simulate_E(Z::Array{Float64, 1})
    @unpack z, T, N, ε, Π_ε, Π_ε_star = Primitives()
    @assert size(Π_ε)[1] == size(Π_ε)[2]

    # X_Z is recession indicator (X_Z[t] == 0 in boom and X_Z[t] == 1 in recession)
    X_Z = Int64.(Z .== z[2])

    # setup the simulation
    X = fill(0, T, N) # allocate memory

    # iterate through each agent
    # X is employment indicator (X[t,i] == 0 if unemployed and X[t,i] == 1 if employed)
    for i in 1:N
        
        # Choose initial employment status based on stationary distribution in initial period.
        X[1, i] = Int64(Π_ε_star[X_Z[1]+ 1]/(Π_ε_star[X_Z[1]+ 1] + Π_ε_star[X_Z[1]+3]) > rand())

        # iterate through future periods
        for t in 2:T
            was_recession = X_Z[t-1]
            was_employed = X[t-1, i]
            is_recession = X_Z[t]

            # Π_ε is constructed as follows:
            # π_{gg11}   π_{gb11}   π_{gg10}   π_{gb10}
            # π_{bg11}   π_{bb11}   π_{bg10}   π_{bb10}
            # π_{gg01}   π_{gb01}   π_{gg00}   π_{gb00}
            # π_{bg01}   π_{bb01}   π_{bg00}   π_{bb00}
            π_e = Π_ε[1 + was_recession + 2*(1-was_employed), 1 + is_recession]
            π_u = Π_ε[1 + was_recession + 2*(1-was_employed), 3 + is_recession]

            X[t, i] = π_e / (π_e + π_u) > rand()

        end
    end
    return X
end

function compute_K(Z::Array{Float64, 1}, a_0::Float64, a_1::Float64, b_0::Float64, b_1::Float64)
    @unpack z = Primitives()

    K = zeros(length(Z))
    K[1] = 11.55 # steady state of complete markets model

    for t = 1:(length(Z) - 1)
        if Z[t] == z[1] 
            K[t+1] = exp(a_0 + a_1 * log(K[t]))
        else
            K[t+1] = exp(b_0 + b_1 * log(K[t]))
        end

        K[t+1] = min(15, max(11, K[t+1]))
    end
    plot(K)
end

####################################################################################################################
################################ Initialize Results structure  #####################################################
####################################################################################################################

function Initialize()
    Z   = simulate_Z()
    E   = simulate_E(Z)
    a_0 = 0.095
    a_1 = 0.965 # chosen so that aggregate capital stays between 11 and 15
    b_0 = 0.085
    b_1 = 0.965 # chosen so that aggregate capital stays between 11 and 15
    # a_0 = 0.0
    # a_1 = 1.0
    # b_0 = 0.0
    # b_1 = 1.0
    K = compute_k(Z, a_0, a_1, b_0, b_1)

    Results(Z, E, a_0, a_1, b_0, b_1, K)
end




