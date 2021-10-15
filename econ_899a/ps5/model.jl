#############################################################################
# ECON 899A Computational Economics
# Problem set 5

# Krusell and Smith (1998)
# Alex von Hafften
# October 13, 2021

# Based on code written by Phil Coyle

# This file contains the functions to estimate the model.
# ./run.jl calls this file and runs the model.
#############################################################################

using Parameters, Interpolations, Optim

####################################################################################################################
##################################### Parameter and result structures ##############################################
####################################################################################################################

# Model primitives
@with_kw struct Primitives
    
    # Preference parameters
    β::Float64                = 0.99                    # discount rate

    # Production technology parameters
    α::Float64                = 0.36                    # elasticity of capital in production
    δ::Float64                = 0.025                   # Capital depreciation
    e_bar::Float64            = 0.3271                  # Labor efficiency per unit of time worked

    # Simulation parameters
    T::Int64                  = 11000                   # Number of periods in simulation
    N::Int64                  = 5000                    # Number of individuals in simulation
    burn_in::Int64            = 1000                    # Number of periods at beginning to throw away before running regression.
end

# Grids for aggregate productivity, employment, aggregate capital, and individual capital
@with_kw struct Grids

    # Aggregate productivity
    z::Array{Float64, 1}      = [1.01, 0.99]            # TFP states
    z_length::Int64           = length(z)               # number of states
    Π_z::Array{Float64}       = get_Π_z()               # aggregate shock transition matrix
    Π_z_star::Array{Float64}  = compute_Π_star(Π_z)     # aggregate shock steady state

    # Employment
    ε::Array{Float64, 1}      = [0, 1]                  # Idiosyncratic employment opportunities
    Π_ε::Array{Float64}       = get_Π_ε()               # employment transition matrix
    Π_ε_star::Array{Float64}  = compute_Π_star(Π_ε)     # stationary distribution of Π

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
    k_grid_srl::Array{Float64, 1} = range(k_min, k_max; length = k_length) # individual capital grid step range length
    k_grid::Array{Float64, 1}     = collect(k_grid_srl)                    # individual capital grid
end

# Results structure
mutable struct Results
    Z::Array{Float64, 1} # Aggregate productivity
    E::Array{Float64, 2} # Employment 
    a_0::Float64         # Intercept in boom capital regression
    a_1::Float64         # Slope in boom capital regression
    b_0::Float64         # Intercept in recession capital regression
    b_1::Float64         # Slope in recession capital regression
    K::Array{Float64, 1} # Aggregate capital
end

# Transition matrices for business cycles and employment.
@with_kw struct Shocks
    # parameters of aggregate state transition matrix:
    d_g::Float64  = 8.0     # Duration (Good Times)
    d_b::Float64  = 8.0     # Duration (Bad Times)

    # transition probabilities for aggregate states
    pgg::Float64 = (d_g-1.0)/d_g
    pgb::Float64 = 1.0 - (d_b-1.0)/d_b
    pbg::Float64 = 1.0 - (d_g-1.0)/d_g
    pbb::Float64 = (d_b-1.0)/d_b

    # Transition matrix for aggregate state
    Π_z::Array{Float64,2} = [pgg pgb; pbg pbb]

    # parameters of aggregate state transition matrix:
    d_ug::Float64 = 1.5     # Unemp Duration (Good Times)
    d_ub::Float64 = 2.5     # Unemp Duration (Bad Times)
    u_g::Float64  = 0.04    # Fraction Unemp (Good Times)
    u_b::Float64  = 0.1     # Fraction Unemp (Bad Times)

    # transition probabilities for aggregate states and staying unemployed
    pgg00::Float64 = (d_ug-1.0)/d_ug
    pbb00::Float64 = (d_ub-1.0)/d_ub
    pbg00::Float64 = 1.25*pbb00
    pgb00::Float64 = 0.75*pgg00

    # transition probabilities for aggregate states and becoming employed
    pgg01::Float64 = (u_g - u_g*pgg00)/(1.0-u_g)
    pbb01::Float64 = (u_b - u_b*pbb00)/(1.0-u_b)
    pbg01::Float64 = (u_b - u_g*pbg00)/(1.0-u_g)
    pgb01::Float64 = (u_g - u_b*pgb00)/(1.0-u_b)

    # transition probabilities for aggregate states and becoming unemployed
    pgg10::Float64 = 1.0 - (d_ug-1.0)/d_ug
    pbb10::Float64 = 1.0 - (d_ub-1.0)/d_ub
    pbg10::Float64 = 1.0 - 1.25*pbb00
    pgb10::Float64 = 1.0 - 0.75*pgg00

    # transition probabilities for aggregate states and staying employed
    pgg11::Float64 = 1.0 - (u_g - u_g*pgg00)/(1.0-u_g)
    pbb11::Float64 = 1.0 - (u_b - u_b*pbb00)/(1.0-u_b)
    pbg11::Float64 = 1.0 - (u_b - u_g*pbg00)/(1.0-u_g)
    pgb11::Float64 = 1.0 - (u_g - u_b*pgb00)/(1.0-u_b)

    # Employment Transition Matrices conditional on aggregate state.
    Π_gg::Array{Float64,2} = [pgg11 pgg01; pgg10 pgg00]
    Π_bg::Array{Float64,2} = [pgb11 pgb01; pgb10 pgb00]
    Π_gb::Array{Float64,2} = [pbg11 pbg01; pbg10 pbg00]
    Π_bb::Array{Float64,2} = [pbb11 pbb01; pbb10 pbb00]

    # Uncoditional employment transition matrix
    Π_ε::Array{Float64,2} = [pgg*Π_gg pgb*Π_gb; pbg*Π_bg pbb*Π_bb]
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

function initialize_k()
    @unpack T = Primitives()

    K = zeros(T)
    K[1] = 11.55
    K
end

####################################################################################################################
################################ Initialize Results structure  #####################################################
####################################################################################################################

function Initialize()
    Z   = simulate_Z()
    E   = simulate_E(Z)
    a_0 = 0.095
    a_1 = 0.999
    b_0 = 0.085
    b_1 = 0.999
    K = zeros()

    Results(Z, E, a_0, a_1, b_0, b_1, K)
end




