# Alex von Hafften
# Problem set 1
# ECON 717B: Applied Econometrics
# Professor Matt Wiswall
# April 7, 2022

# This file holds functions that are called in ./run.jl

using Parameters, Distributions, Random, Optim

# parameters structure
mutable struct parameters
    π_1::Float64
    π_2::Float64
    μ_1::Float64
    μ_2::Float64
    σ_1::Float64
    σ_2::Float64
    ρ::Float64
end

# data object
mutable struct data
    N::Int64             # number of trials
    ε_1::Vector{Float64} # idiosyncratic component for skill in occupation 1
    ε_2::Vector{Float64} # idiosyncratic component for skill in occupation 2
    S_1::Vector{Float64} # skill in occupation 1
    S_2::Vector{Float64} # skill in occupation 2
    W_1::Vector{Float64} # potential wages in occupation 1
    W_2::Vector{Float64} # potential wages in occupation 2
    D::Vector{Int64}     # occupation choice (observed)
    W::Vector{Float64}   # wages (observed)
end

# simulates data based on the parameters and seed
function simulate_data(θ::parameters; seed = -1)

    # trials
    N = 1000

    # covariance matrix
    σ_11 = θ.σ_1 * θ.σ_1
    σ_22 = θ.σ_2 * θ.σ_2
    σ_12 = θ.ρ * θ.σ_1 * θ.σ_2

    # create ε-distribution
    ε_distribution = MvNormal([θ.μ_1, θ.μ_2], [σ_11 σ_12; σ_12 σ_22])

    # set seed if it is specified
    if seed != -1
        Random.seed!(seed)
    end

    # draw ε_k
    ε = rand(ε_distribution, N)
    ε_1 = ε[1,:]
    ε_2 = ε[2,:]

    # compute skills
    S_1 = exp.(θ.π_1 .+ ε_1)
    S_2 = exp.(θ.π_2 .+ ε_2)

    # compute wages
    W_1 = θ.π_1 .* S_1
    W_2 = θ.π_2 .* S_2

    # occupation choice
    D = W_1 .>= W_2

    # observed wages
    W = D .* W_1 .+ (1 .- D) .* W_2

    data(N, ε_1, ε_2, S_1, S_2, W_1, W_2, D, W)
end

# observation level likelihood
function likelihood(θ::parameters, D::Int64, W::Float64)

    # if observation chooses occupation 1
    if D == 1

        # i denotes parameter of chosen occupation
        # j denotes parameter of other occupation

        π_i = θ.π_1
        π_j = θ.π_2
        μ_i = θ.μ_1
        μ_j = θ.μ_2
        σ_i = θ.σ_1
        σ_j = θ.σ_2
        ρ   = θ.ρ

    # if observation chooses occupation 2
    elseif D == 0

        # i denotes parameter of chosen occupation
        # j denotes parameter of other occupation

        π_i = θ.π_2
        π_j = θ.π_1
        μ_i = θ.μ_2
        μ_j = θ.μ_1
        σ_i = θ.σ_2
        σ_j = θ.σ_1
        ρ   = θ.ρ

    else
        error("D is invalid value")
    end

    ε_i = log(W_i/π_i) - μ_i

    cdf(Normal(μ_j + (σ_j/σ_i) * ρ * (ε_i - μ_i), (1-ρ^2) * σ_j^2), log(W_i) - μ_j)
end

# log-likelihood
function log_likelihood(θ::parameters, D::data)

    result = 0.0
    for i = 1:D.N
        result += log(likelihood(θ, D.D[i], D.W[i]))
    end

    result
end

# estimate μ_1 and ρ using maximum likelihood
function mle(θ_0::parameters, D::data)

    function obj(guess)

        ρ_hat = 2 * (exp(guess[2])/(1 + exp(guess[2]))) - 1

        θ_hat = parameters(θ_0.π_1, θ_0.π_2, guess[1], θ_0.μ_2, θ_0.σ_1, θ_0.σ_2, ρ_hat)

        -log_likelihood(θ_hat, D)
    end

    opt = optimize(obj, [θ_0.μ_1, θ_0.ρ]).minimizer

    return(opt[1], 2 * (exp(opt[2])/(1 + exp(opt[2]))) - 1)
end

# identification plot for μ_1
function id_plot_μ_1(μ_1_hat::Float64, ρ_hat::Float64, θ_0::parameters, D::data)

    μ_1 = 0.0:0.1:2.0

    ll_μ_1_ρ_0 = zeros(length(μ_1))
    ll_μ_1_ρ_hat = zeros(length(μ_1))

    for i in 1:length(μ_1)
        ll_μ_1_ρ_0[i] =  log_likelihood(parameters(θ_0.π_1, θ_0.π_2, μ_1[i], θ_0.μ_2, θ_0.σ_1, θ_0.σ_2, θ_0.ρ), D)
        ll_μ_1_ρ_hat[i] =  log_likelihood(parameters(θ_0.π_1, θ_0.π_2, μ_1[i], θ_0.μ_2, θ_0.σ_1, θ_0.σ_2, ρ_hat), D)
    end

    plot(μ_1, ll_μ_1_ρ_0, label = "Log Likelihoof for ρ = ρ_0")
    plot!(μ_1, ll_μ_1_ρ_hat, label = "Log Likelihoof for ρ = ρ_hat")
    vline!([μ_1_hat], label = "μ_1_hat")
    vline!([θ_0.μ_1], label = "μ_1_0", legend = :bottomleft)
    xlabel!("μ_1")
    ylabel!("Log Likelihood")
end


# identification plot for ρ
function id_plot_ρ(μ_1_hat::Float64, ρ_hat::Float64, θ_0::parameters, D::data)

    ρ = -1.0:0.1:1.0

    ll_ρ_μ_1_0 = zeros(length(ρ))
    ll_ρ_μ_1_hat = zeros(length(ρ))

    for i in 1:length(ρ)
        ll_ρ_μ_1_0[i] =  log_likelihood(parameters(θ_0.π_1, θ_0.π_2, θ_0.μ_1, θ_0.μ_2, θ_0.σ_1, θ_0.σ_2, ρ[i]), D)
        ll_ρ_μ_1_hat[i] =  log_likelihood(parameters(θ_0.π_1, θ_0.π_2, μ_1_hat, θ_0.μ_2, θ_0.σ_1, θ_0.σ_2, ρ[i]), D)
    end

    plot(ρ, ll_ρ_μ_1_0, label = "Log-likelihood for μ_1 = μ_1_0")
    plot!(ρ, ll_ρ_μ_1_hat, label = "Log-likelihood for μ_1 = μ_1_hat")
    vline!([ρ_hat], label = "ρ_hat")
    vline!([θ_0.ρ], label = "ρ_0", legend = :bottomleft)
    xlabel!("ρ")
    ylabel!("Log Likelihood")
end