####################################################################################################################
# ECON 899A Computational Economics
# Problem set 5

# Krusell and Smith (1998)
# Alex von Hafften
# October 13, 2021

# Based on code written by Phil Coyle

# This file initializes the results structure with business cycles and employment.
####################################################################################################################

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
    a0::Float64         # Intercept in boom capital regression
    a1::Float64         # Slope in boom capital regression
    b0::Float64         # Intercept in recession capital regression
    b1::Float64         # Slope in recession capital regression

    # Other
    K::Array{Float64, 1}  # Aggregate capital
    R2::Float64           # Explanatory power of forecasting regression
end

####################################################################################################################
####### Functions to simulate business cycles and employment based on business cycles. #############################
####################################################################################################################

# simulate business cycles
function simulate_Z(;seed = missing)
    @unpack T = Simulation()
    @unpack Π_z = Shocks()

    # setup the simulation
    Z = fill(0, T) # allocate memory
    Z[1] = 1

    # Optional seed
    if (!ismissing(seed))
        Random.seed!(seed)
    end
    dist = Uniform(0, 1)
   
    # iterate through each period
    # X is recession indicator (X[t] == 0 in boom and X[t] == 1 in recession)
    for t in 2:T  
        shock = rand(dist)

        if Z[t-1] == 1 && shock < Π_z[1, 1]
            Z[t] = 1
        elseif Z[t-1] == 1 && shock > Π_z[1, 1]
            Z[t] = 2
        elseif Z[t-1] == 2 && shock > Π_z[2, 2]
            Z[t] = 1
        else # elseif Z[t-1] == 2 && shock < Π_z[2, 2]
            Z[t] = 2
        end
    end
    
    return Z
end

# simulate employment based on business cycles
function simulate_E(Z::Array{Int64, 1}; seed = missing)
    @unpack T, N = Simulation()
    @unpack Π_gg, Π_gb, Π_bg, Π_bb = Shocks()

    # setup the simulation
    E = fill(0, N, T) # allocate memory
    E[:,1] .= 1

    # Optional seed
    if (!ismissing(seed))
        Random.seed!(seed)
    end
    dist = Uniform(0, 1)

    # iterate through each agent and each period
    for i in 1:N
        for t in 2:T
            
            # Choose employment transition probabilities based on aggregate states
            if Z[t-1] == 1 && Z[t] == 1
                p11 = Π_gg[1,1]
                p00 = Π_gg[2,2]
            elseif Z[t-1] == 1 && Z[t] == 2
                p11 = Π_gb[1,1]
                p00 = Π_gb[2,2]
            elseif Z[t-1] == 2 && Z[t] == 1
                p11 = Π_bg[1,1]
                p00 = Π_bg[2,2]
            else # elseif Z[t-1] == 2 && Z[t] == 2
                p11 = Π_bb[1,1]
                p00 = Π_bb[2,2]
            end

            # Draw shock
            shock = rand(dist)

            # Assign current employment based on previous employment
            if E[i,t-1] == 1 && shock < p11
                E[i,t] = 1
            elseif E[i,t-1] == 1 && shock > p11
                E[i,t] = 2
            elseif E[i,t-1] == 2 && shock < p00
                E[i,t] = 2
            else # elseif E[i,t-1] == 2 && shock > p00
                E[i,t] = 1
            end
        end
    end
    return E
end

####################################################################################################################
################################ Initialize Results structure  #####################################################
####################################################################################################################

function Initialize()
    @unpack max_iterations, T = Simulation()
    @unpack n_k, n_K, n_ε, n_z, k_grid = Grids()

    Z   = simulate_Z(seed = 12003030)
    E   = simulate_E(Z; seed = 12003030)

    value_function  = zeros(n_k, n_ε, n_K, n_z)
    policy_function = zeros(n_k, n_ε, n_K, n_z)
    
    # guess based on previous model runs
    a0 = 0.11914580022527435
    a1 = 0.9471228553642084
    b0 = 0.2029256394484849
    b1 = 0.9133652265162197
    
    K   = zeros(T)
    R2  = 0

    Results(Z, E, value_function, policy_function, a0, a1, b0, b1, K, R2)
end
