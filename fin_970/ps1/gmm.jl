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
    beta_1::Array{Float64}     # point estimates
    S_1::Array{Float64}     # newey-West variance-covariance matrix esimate
    se_1::Array{Float64}    # newey-West standard errors

    # second stage
    W::Array{Float64}       # optimal weighting matrix
    beta_2::Array{Float64}     # point estimates
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
    function gmm_obj(beta::Array{Float64}, W)
        moments = 1/T*Z'*(Y .- X*beta)
        return moments'*W*moments
    end

    # first stage coefficient estimates
    beta_1 = optimize(beta -> gmm_obj(beta, I), zeros(k)).minimizer
    
    # newey-west variance covariance estimate
    function newey_west_S(beta::Array{Float64}, L::Int64)
        
        h = (Y.-X*beta).*Z
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
    S_1 = newey_west_S(beta_1, L)
    se_1 = sqrt.(diag(S_1))

    # optimal weighting matrix
    W = inv(S_1)

    # second stage point estimates
    beta_2 = optimize(beta -> gmm_obj(beta, W), zeros(k)).minimizer

    # second stage newey west var-cov matrix and ses.
    S_2 = newey_west_S(beta_2, L)
    se_2 = sqrt.(diag(S_2))

    # R-squared
    R_2 = 1 - var(Y .- X*beta_2)/var(Y)

    return GMM(Y, X, Z, T, r, k, L, beta_1, S_1, se_1, W, beta_2, S_2, se_2, R_2)
end

