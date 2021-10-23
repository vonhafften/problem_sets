##########################################################################################
# Problem set 7
# ECON 899A: Computational Economics
#
# Alex von Hafften
# October 27, 2021
#
# This code is called by run.jl
##########################################################################################

using Parameters, Random, Distributions, Plots, LinearAlgebra, Interpolations, Optim

##########################################################################################
#### Utility functions
##########################################################################################

# function to draw true shocks
function shocks(T::Int64, H::Int64; seed = missing)
    if !ismissing(seed)
        Random.seed!(seed)
    end
    rand(Normal(0, 1), T, H) # Draw shocks
end

# computes x based on shocks,  ρ, and σ
function dgp(ε::Array{Float64, 2}, ρ::Float64, σ::Float64)
    
    T = size(ε)[1] # get dimensions
    H = size(ε)[2]

    x = zeros(T, H) # initialize x matrix

    # fill in x based on parameters
    for i_H = 1:H
        x[1, i_H] = ε[1, i_H]
        for t = 2:T
            x[t, i_H] = ρ * x[t - 1, i_H] + σ * ε[t, i_H]
        end
    end
    x
end

# transforms x in 3d array
# first layer of third dimenion is just x (for sample moment = mean)
# second layer of third dimension is (x - x_bar)^2 (for sample moment = variance)
# third layer of third dimension is (x - x_bar)*(x_lag - x_bar) (for sample moment = first order autocorrelation)
function m(x::Array{Float64, 2})

    T = size(x)[1] # get dimensions
    H = size(x)[2]

    x_bar = repeat(sum(x; dims = 1)./T, T) # compute mean by simulation
    x_lag = vcat(zeros(H)', x[1:T-1, :])   # lag x

    m = zeros(T, H, 3)
    m[:,:,1] = x
    m[:,:,2] = (x .- x_bar).^2
    m[:,:,3] = (x .- x_bar) .* (x_lag .- x_bar)

    m
end

# computes sample average of transformed x
function M(m::Array{Float64, 3})
    T = size(m)[1] # get dimensions
    H = size(m)[2]
    
    # computes sample average
    M = reshape(sum(m; dims = 1:2)/(T*H), 3)

end

# Create Primitive structure to hold true data and other parameters
@with_kw struct Primitives

    # True data
    ρ_0::Float64           = 0.5                       # True ρ parameter
    σ_0::Float64           = 1.0                       # True σ paramter
    T::Int64               = 200                       # Length of data/simulations
    ε::Array{Float64, 2}   = shocks(T, 1; seed = 1234) # true data error
    x_0::Array{Float64, 2} = dgp(ε, ρ_0, σ_0)          # true data
    m_0::Array{Float64, 3} = m(x_0)                    # true transformed data
    M_T::Array{Float64, 1} = M(m_0)                    # true data moments

    # Simulation parameters
    H::Int64               = 10  # Number of simulations
    i_T::Int64             = 4   # Number of lags

    # Grids for plotting
    ρ_grid::Array{Float64, 1} = 0.35:0.01:0.65
    σ_grid::Array{Float64, 1} = 0.8:0.01:1.2
end

# results structure with which moments to consider, shocks, first state estimates,
# estimate of optimal weighting matrix, and second stage estimates.
mutable struct Results
    which_moments::Array{Int64, 1}  # Which moments to include in objective function
    e::Array{Float64, 2}            # Simulated shocks
    ρ_hat_1::Float64                # First stage ρ estimator
    σ_hat_1::Float64                # First stage σ estimator 
    W_hat::Array{Float64, 2}   # Estimator of optimal weight matrix
    ρ_hat_2::Float64                # First stage ρ estimator
    σ_hat_2::Float64                # First stage σ estimator 
end

# Initialize results structure
function Initialize(which_moments; seed)
    @unpack T, H = Primitives()

    e = shocks(T, H)
    ρ_hat_1 = 0
    σ_hat_1 = 0
    W_hat = zeros(length(which_moments), length(which_moments))
    ρ_hat_2 = 0
    σ_hat_2 = 0

    Results(which_moments, e, ρ_hat_1, σ_hat_1, W_hat, ρ_hat_2, σ_hat_2)
end

# Calculate objective function for parameters b, 
# shocks e, moments which_moments, weighting W, and population data moments M_T
function J_TH(b, e, which_moments, W, M_T)
    ρ = b[1] # unpack b
    σ = b[2]
    
    M_TH = M(m(dgp(e, ρ, σ))) # compute M_TH

    # evaluate J_TH
    (M_T[which_moments] - M_TH[which_moments])' * W * (M_T[which_moments] - M_TH[which_moments])
end

# Plot objective function surface
function plot_J_TH_surface(R::Results, stage::Int64)
    @unpack ρ_grid, σ_grid, M_T = Primitives()

    # choose weighting matrix based on stage argument
    if stage == 1
        W = I
    elseif stage == 2
        W = R.W_hat
    else 
        error("Specify valid stage argument: 1 or 2.")
    end

    # create plotting function
    function J_TH_plotting(ρ, σ)
        return J_TH([ρ, σ], R.e, R.which_moments, W, M_T)
    end

    # Surface of J_TH plot
    surface(ρ_grid, σ_grid, J_TH_plotting);
    xlabel!("ρ");
    ylabel!("σ")
end

# estimate coefficients 
function b_hat(R::Results, stage::Int64)
    @unpack ρ_0, σ_0, M_T = Primitives()

    # choose weighting matrix based on stage argument
    if stage == 1
        W = I
    elseif stage == 2
        W = R.W_hat
    else 
        error("Specify valid stage argument: 1 or 2.")
    end

    opt = optimize(b -> J_TH(b, R.e, R.which_moments, W, M_T), [ρ_0, σ_0])
    Optim.minimizer(opt)
end

# estimate optimal weighting matrix
function W_hat(R::Results)
    @unpack T, H, i_T = Primitives()

    # compute M_TH
    m_1 = m(dgp(R.e, R.ρ_hat_1, R.σ_hat_1))
    M_TH_1 = M(m_1)

    function compute_Γ(j)
        Γ = zeros(3, 3)

        for i_1 = 1:3
            for i_2 in 1:3
                Γ[i_1, i_2] = sum((m_1[j+1:T,:,i_1] .- M_TH_1[i_1]) * (m_1[1:T-j, :, i_2] .- M_TH_1[i_2])')/ (T * H)
            end
        end
        Γ[R.which_moments, R.which_moments]
    end

    S_TH = compute_Γ(0)

    for j = 1:i_T
        Γ_j = compute_Γ(j)
        S_TH += (1 - (j/(i_T + 1))) * (Γ_j + Γ_j')
    end

    S_TH = (1 + 1/H) * S_TH

    inv(S_TH)
end

function estimate(which_moments; seed = missing)
    # Initialize results object
    R = Initialize(which_moments; seed)

    # First stage
    b_hat_1 = b_hat(R, 1)
    R.ρ_hat_1 = b_hat_1[1]
    R.σ_hat_1 = b_hat_1[2]

    # Estimate optimal weighting matrix
    R.W_hat = W_hat(R)

    # second stage 
    b_hat_2 = b_hat(R, 2)
    R.ρ_hat_2 = b_hat_2[1]
    R.σ_hat_2 = b_hat_2[2]

    R
end