#############################################################################
# ECON 899A Computational Economics
# Problem set 5

# Krusell and Smith (1998)
# Alex von Hafften
# October 13, 2021

# Based on code written by Phil Coyle

# This file initializes the results structure with business cycles and employment.
#############################################################################

include("01_parameters.jl");

# Results structure
mutable struct Results

    # Realization of business cycles and employment
    Z::Array{Int64, 1} # Aggregate productivity state
    E::Array{Int64, 2} # Employment 

    # Value and policy functions
    value_function::Array{Float64,4}
    policy_function::Array{Float64,4}

    # Coefficients on forecasting regression
    a_0::Float64         # Intercept in boom capital regression
    a_1::Float64         # Slope in boom capital regression
    b_0::Float64         # Intercept in recession capital regression
    b_1::Float64         # Slope in recession capital regression

    # Other
    K::Array{Float64, 1} # Aggregate capital
    R2::Float64          # Explanatory power of forecasting regression
end

####################################################################################################################
####### Functions to simulate business cycles and employment based on business cycles. #############################
####################################################################################################################

# simulate business cycles
function simulate_Z()
    @unpack T = Simulation()
    @unpack Π_z = Shocks()
    @assert size(Π_z)[1] == size(Π_z)[2] # square required

    # setup the simulation
    Z = fill(0, T) # allocate memory
   
    # iterate through each period
    # X is recession indicator (X[t] == 0 in boom and X[t] == 1 in recession)
    for t in 2:T  
        Z[t] = Π_z[Z[t-1]+1, 1] < rand()
    end
    
    return Z
end

# simulate employment based on business cycles
function simulate_E(Z::Array{Float64, 1})
    @unpack T, N = Simulation()
    @unpack Π_ε = Shocks()
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





