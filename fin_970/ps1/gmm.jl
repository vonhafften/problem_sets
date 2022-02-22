# Alex von Hafften 
# FIN 970: Asset Pricing
# HW 1 Problem 1
# Professor Ivan Shaliastovich

# This code estimates a linear regression via gmm.
# Produces point estimates, Newey-West standard errors, and R^2

using Optim, LinearAlgebra

# struct to hold gmm estimation
mutable struct GMM
    # data 
    Y::Array{Float64}       # dependent variable
    X::Array{Float64}       # independent variables
    Z::Array{Float64}       # instruments

    # sizes
    T::Int64                 # number of observations
    r::Int64                 # number of moments
    k::Int64                 # number of parameters
    L::Int64                 # number of lags for newey-west SEs

    # first stage
    β_1::Array{Float64}     # point estimates
    S_1::Array{Float64}     # newey-West variance-covariance matrix esimate
    se_1::Array{Float64}    # newey-West standard errors

    # second stage
    W::Array{Float64}       # optimal weighting matrix
    β_2::Array{Float64}     # point estimates
    S_2::Array{Float64}     # newey-West variance-covariance matrix esimate
    se_2::Array{Float64}    # newey-West standard errors

    # post estimation
    R_2::Float64             # r-squared
end

# estimate linear regression with gmm
function gmm(Y::Array{Float64}, X::Array{Float64}, Z::Array{Float64}, L::Int64)
    
    # get parameters
    T = size(X)[1]
    k = size(X)[2]
    r = size(Z)[2]

    # gmm objective function
    function gmm_obj(β::Array{Float64}, W)
        moments = 1/T*Z'*(Y .- X*β)
        return moments'*W*moments
    end

    # first stage coefficient estimates
    β_1 = optimize(β -> gmm_obj(β, I), zeros(k)).minimizer
    
    # newey-west variance covariance estimate
    function newey_west_S(β::Array{Float64}, L::Int64)
        
        h = (Y.-X*β).*Z
        S = 1/T * h' * h

        for i = 1:L
            h_1 = h[(i+1):end, :]
            h_2 = h[1:(end-i), :]
            G_q = 1/T .* h_1' * h_2
            S += (1 - i/(L+1)) .* (G_q + G_q')
        end

        return S
    end

    # first stage newey west var-cov matrix and ses.
    S_1 = newey_west_S(β_1, L)
    se_1 = sqrt.(diag(S_1))

    # optimal weighting matrix
    W = inv(S_1)

    # second stage point estimates
    β_2 = optimize(β -> gmm_obj(β, W), zeros(k)).minimizer

    # second stage newey west var-cov matrix and ses.
    S_2 = newey_west_S(β_2, L)
    se_2 = sqrt.(diag(S_2))

    # R-squared
    R_2 = 1 - var(Y .- X*β_2)/var(Y)

    return GMM(Y, X, Z, T, r, k, L, β_1, S_1, se_1, W, β_2, S_2, se_2, R_2)
end

