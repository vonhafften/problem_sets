using Plots, Tables, DataFrames, CSV, StatFiles, Optim

df = DataFrame(load("Mortgage_performance_data.dta"))

# create new columns
df[!, :i_open_year2_i_open_year5] = df[!, :i_open_year2] .+ df[!, :i_open_year5]
df[!, :c] .= 1

# Define y and x
df_y = select(df, :i_close_first_year)
df_x = select(df, :c, :i_large_loan, :i_medium_loan, :rate_spread,
              :i_refinance, :age_r, :cltv, :dti, :cu,  :first_mort_r,
              :score_0, :score_1, :i_FHA, :i_open_year2_i_open_year5)

y = Matrix(df_y)
x = Matrix(df_x)

N = size(x)[1]
K = size(x)[2]

# Log-likelihood
function log_likelihood(β::Array{Float64, 1})
    L = 0
    for i = 1:N 
        Λ_i = exp(β' * x[i,:])/(1 + exp(β' * x[i,:]))
        L += Λ_i^y[i] * (1 - Λ_i)^(1 - y[i])
    end
    return L
end

# negative of log-likelihood for using builtin optimizing packages
function log_likelihood_neg(β::Array{Float64, 1})
    -log_likelihood(β)
end

# Score
function score(β::Array{Float64, 1})
    g = zeros(K)
    for i = 1:N 
        Λ_i = exp(β' * x[i,:])/(1 + exp(β' * x[i,:]))
        g += (y[i] - Λ_i).*x[i,:]
    end
    return g  
end

# negative of score for using builtin optimizing packages
function score_neg!(storage::Array{Float64, 1}, β::Array{Float64, 1})
    storage .= -score(β)
end

# Hessian
function hessian(β::Array{Float64, 1})
    H = zeros(K, K)
    for i = 1:N 
        Λ_i = exp(β' * x[i,:])/(1 + exp(β' * x[i,:]))
        H .+= -Λ_i * (1 - Λ_i) * x[i,:] * x[i,:]'
    end
    return H   
end

# negative of hessian for using builtin optimizing packages
function hessian_neg!(storage::Array{Float64, 1}, β::Array{Float64, 1})
    storage .= -hessian(β)
end

# numerical calculate the score
function score_numerical(β::Array{Float64, 1})
    ε = 1e-10
    g = zeros(K)
    for i = 1:K
        # adding ε to β_i
        β_upper = copy(β)
        β_upper[i] += ε
        L_upper = log_likelihood(β_upper)

        # subtracting ε to β_i
        β_lower = copy(β)
        β_lower[i] -= ε
        L_lower = log_likelihood(β_lower)

        # central approximation
        g[i] = (L_upper - L_lower)/(2*ε)
    end
    return g
end

function hessian_numerical(β::Array{Float64, 1})
    ε = 1e-10
    H = zeros(K, K)
    for i = 1:K
        # adding ε to β_i
        β_upper = copy(β)
        β_upper[i] += ε
        g_upper = score_numerical(β_upper)

        # subtracting ε to β_i
        β_lower = copy(β)
        β_lower[i] -= ε
        g_lower = score_numerical(β_lower)

        # central approximation
        H[i, :] = (g_upper .- g_lower)./(2*ε)
    end
    return H
end


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

# Exercise 1
β_1 = prepend!(zeros(K-1), -1)
log_likelihood(β_1)
score(β_1)
hessian(β_1)

# Exercise 2
score_numerical(β_1)
hessian_numerical(β_1)

# Exercise 3
newton(log_likelihood, score, hessian, β_1)

# Exercise 4
minimization_nelder_mead = optimize(log_likelihood_neg, β_1, NelderMead())
minimization_lbfgs = optimize(log_likelihood_neg, score_neg!, β_1, LBFGS())
minimization_newton = optimize(log_likelihood_neg, score_neg!, hessian_neg!, β_1, Newton())