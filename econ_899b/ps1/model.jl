# Computational Economics
# Professor JF Houde
# Problem set 1
# Alex von Hafften 
# November 7, 2021

# Log-likelihood
function log_likelihood(β::Array{Float64, 1}, x::Matrix{Float64}, y::Matrix{Float64})
    N = size(x)[1]
    L = 0
    for i = 1:N 
        Λ = exp(β' * x[i,:])/(1 + exp(β' * x[i,:]))
        L += log((Λ^y[i]) * ((1 - Λ)^(1 - y[i])))
    end
    return L
end

# Score
function score(β::Array{Float64, 1}, x::Matrix{Float64}, y::Matrix{Float64})
    N = size(x)[1]
    g = zeros(K)
    for i = 1:N 
        Λ = exp(β' * x[i,:])/(1 + exp(β' * x[i,:]))
        g += (y[i] - Λ).*x[i,:]
    end
    return g  
end

# Hessian
function hessian(β::Array{Float64, 1}, x::Matrix{Float64}, y::Matrix{Float64})
    N = size(x)[1]
    K = size(x)[2]

    H = zeros(K, K)
    for i = 1:N 
        Λ = exp(β' * x[i,:])/(1 + exp(β' * x[i,:]))
        H .+= -Λ * (1 - Λ) * x[i,:] * x[i,:]'
    end
    return H   
end

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

# numerical calculate the score
function score_numerical(β::Array{Float64, 1}, x::Matrix{Float64}, y::Matrix{Float64})
    N = size(x)[1]
    
    ε = 1e-10
    g = zeros(K)
    for i = 1:K
        # adding ε to β_i
        β_upper = copy(β)
        β_upper[i] += ε
        L_upper = log_likelihood(β_upper, x, y)

        # subtracting ε to β_i
        β_lower = copy(β)
        β_lower[i] -= ε
        L_lower = log_likelihood(β_lower, x, y)

        # central approximation
        g[i] = (L_upper - L_lower)/(2*ε)
    end
    return g
end

function hessian_numerical(β::Array{Float64, 1}, x::Matrix{Float64}, y::Matrix{Float64})
    ε = 1e-10
    H = zeros(K, K)
    for i = 1:K
        # adding ε to β_i
        β_upper = copy(β)
        β_upper[i] += ε
        g_upper = score_numerical(β_upper, x, y)

        # subtracting ε to β_i
        β_lower = copy(β)
        β_lower[i] -= ε
        g_lower = score_numerical(β_lower, x, y)

        # central approximation
        H[i, :] = (g_upper .- g_lower)./(2*ε)
    end
    return H
end

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

# Newton's method
# f is objective function
# g is first derivative (gradient)
# h is second derivative (Hessian)
# z... are other arguments for f, g, and h
function newton(f, g, h, guess, z...)
    γ = 1 # adjustment parameter

    err, i = 100, 0

    while err > 10e-12
        i += 1

        dx = g(guess, z...)
        ddx =  h(guess, z...)

        next_guess = guess .- inv(ddx) * dx
        
        # apply sup norm for convergence
        err = maximum(abs.(guess - next_guess))
        guess = (1-γ) .* guess .+ γ .* next_guess
    end

    guess, f(guess, z...), i
end

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

# negative of score for using builtin optimizing packages
function score!(storage::Array{Float64, 1}, β::Array{Float64, 1})
    storage .= (-1) .* score(β, x, y)
end

# negative of hessian for using builtin optimizing packages
function hessian!(storage::Array{Float64, 2}, β::Array{Float64, 1})
    storage .= (-1) .* hessian(β, x, y)
end