# Alex von Hafften
# Problem set 1
# ECON 717B: Applied Econometrics
# Professor Matt Wiswall
# April 7, 2022

# This file holds functions that are called in ./run.jl

using Parameters, Distributions, Random, Optim, Statistics, KernelDensity, DataFrames

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
function simulate_data(θ::parameters; seed = -1, N = 1000)

    # covariance matrix
    σ_11 = θ.σ_1 * θ.σ_1
    σ_22 = θ.σ_2 * θ.σ_2
    σ_12 = θ.ρ * θ.σ_1 * θ.σ_2

    # create ε-distribution
    ε_distribution = MvNormal([0, 0], [σ_11 σ_12; σ_12 σ_22])

    # set seed if it is specified
    if seed != -1
        Random.seed!(seed)
    end

    # draw ε_k
    ε = rand(ε_distribution, N)
    ε_1 = ε[1,:]
    ε_2 = ε[2,:]

    # compute skills
    S_1 = exp.(θ.μ_1 .+ ε_1)
    S_2 = exp.(θ.μ_2 .+ ε_2)

    # compute wages
    W_1 = θ.π_1 .* S_1
    W_2 = θ.π_2 .* S_2

    # occupation choice
    D = W_1 .>= W_2

    # observed wages
    W = D .* W_1 .+ (1 .- D) .* W_2

    data(N, ε_1, ε_2, S_1, S_2, W_1, W_2, D, W)
end



# computes observation-level likelihood
# using built-in Normal cdf approximation
function likelihood(θ::parameters, D::Int64, W::Float64)
    if D == 1
        return pdf(Normal(), (log(W/θ.π_1) - θ.μ_1)/θ.σ_1) *
            cdf(Normal(), (log(θ.π_1) - log(θ.π_2) +θ.μ_1 - θ.μ_2 +(1-(θ.σ_2/θ.σ_1*θ.ρ))*(log(W/θ.π_1)-θ.μ_1))/(sqrt((1-θ.ρ^2))*θ.σ_2))
    elseif D == 0
        return pdf(Normal(), (log(W/θ.π_2) - θ.μ_2)/θ.σ_2) *
            cdf(Normal(), (log(θ.π_2) - log(θ.π_1) +θ.μ_2 - θ.μ_1 +(1-(θ.σ_1/θ.σ_2*θ.ρ))*(log(W/θ.π_2)-θ.μ_2))/(sqrt((1-θ.ρ^2))*θ.σ_1))
    end
end

# computes log-likelihood for data
function log_likelihood(θ::parameters, d::data; method = "builtin")

    result = 0.0

    # uses built-in normal cdf approximation
    if method == "builtin"
        # iterates over observations
        for i = 1:d.N
            result += log(likelihood(θ, d.D[i], d.W[i]))
        end
    
    # uses simulation to evaluate likelihood
    elseif method == "simulation"
        
        # create simulated data
        simulations = simulate_data(θ; seed = 123, N = 10000)
        
        # estimate kernel density for wage distribution for each occupation
        w_kde_1 = InterpKDE(kde(simulations.W[Bool.(simulations.D .== 1)]))
        w_kde_2 = InterpKDE(kde(simulations.W[Bool.(simulations.D .== 0)]))

        # get fraction in each occupation
        N_1 =  sum((simulations.D .== 1)) / simulations.N
        N_2 =  sum((simulations.D .== 0)) / simulations.N

        # iterates over observations
        for i = 1:d.N
            if d.D[i] == 1
                result += log(pdf(w_kde_1, d.W[i]) * N_1)
            elseif d.D[i] == 0
                result += log(pdf(w_kde_2, d.W[i]) * N_2)
            end
        end
    end

    result
end

# estimates μ_1 and ρ using maximum likelihood
function mle(θ_0::parameters, d::data; method = "builtin")

    # the objective function to minimizer for the optimizer
    function obj(guess)

        # ρ is confined to be between -1 and 1.
        # here, I use a rescaled logistic function, so that the optimizer can search with bounds.
        ρ_hat = 2 * (exp(guess[2])/(1 + exp(guess[2]))) - 1

        θ_hat = parameters(θ_0.π_1, θ_0.π_2, guess[1], θ_0.μ_2, θ_0.σ_1, θ_0.σ_2, ρ_hat)

        -log_likelihood(θ_hat, d; method = method)
    end

    opt = optimize(obj, [θ_0.μ_1, θ_0.ρ]).minimizer

    return(opt[1], 2 * (exp(opt[2])/(1 + exp(opt[2]))) - 1)
end

# 2d identification plot for μ_1 and ρ
function id_plot_3d(μ_1_hat::Float64, ρ_hat::Float64, θ_0::parameters, d::data; method = "builtin")

    # grids
    μ_1 = 1.0:0.05:1.5
    ρ = 0.0:0.05:0.5

    # function for plotting
    plot_f(μ_1_i, ρ_i) =  log_likelihood(parameters(θ_0.π_1, θ_0.π_2, μ_1_i, θ_0.μ_2, θ_0.σ_1, θ_0.σ_2, ρ_i), d; method=method)

    # create plot
    plot(μ_1, ρ, plot_f, st = :contour, fill=true)
    scatter!([θ_0.μ_1], [θ_0.ρ], label = "(μ_1_0, ρ_0)")
    scatter!([μ_1_hat], [ρ_hat], label = "(μ_1_hat, ρ_hat)")
    xlabel!("μ_1")
    ylabel!("ρ")
end

# identification plot for μ_1
function id_plot_μ_1(μ_1_hat::Float64, ρ_hat::Float64, θ_0::parameters, d::data; method = "builtin")

    # μ_1 grid
    μ_1 = 1.0:0.01:1.5

    # initialize object for log-likelihood vector
    ll_μ_1_ρ_0 = zeros(length(μ_1))
    ll_μ_1_ρ_hat = zeros(length(μ_1))

    # evaluate log-likelihood across grid
    for i in 1:length(μ_1)
        ll_μ_1_ρ_0[i] =  log_likelihood(parameters(θ_0.π_1, θ_0.π_2, μ_1[i], θ_0.μ_2, θ_0.σ_1, θ_0.σ_2, θ_0.ρ), d; method=method)
        ll_μ_1_ρ_hat[i] =  log_likelihood(parameters(θ_0.π_1, θ_0.π_2, μ_1[i], θ_0.μ_2, θ_0.σ_1, θ_0.σ_2, ρ_hat), d; method=method)
    end

    # plot log-likelihood
    plot(μ_1, ll_μ_1_ρ_0, label = "Log Likelihood for ρ = ρ_0")
    plot!(μ_1, ll_μ_1_ρ_hat, label = "Log Likelihood for ρ = ρ_hat")
    vline!([μ_1_hat], label = "μ_1_hat")
    vline!([θ_0.μ_1], label = "μ_1_0", legend = :bottomleft)
    xlabel!("μ_1")
    ylabel!("Log Likelihood")
end

# identification plot for ρ
function id_plot_ρ(μ_1_hat::Float64, ρ_hat::Float64, θ_0::parameters, d::data; method = "builtin")

    # ρ grid
    ρ = 0.0:0.01:0.5

    # initialize object for log-likelihood vector
    ll_ρ_μ_1_0 = zeros(length(ρ))
    ll_ρ_μ_1_hat = zeros(length(ρ))

    # evaluate log-likelihood across grid
    for i in 1:length(ρ)
        ll_ρ_μ_1_0[i] =  log_likelihood(parameters(θ_0.π_1, θ_0.π_2, θ_0.μ_1, θ_0.μ_2, θ_0.σ_1, θ_0.σ_2, ρ[i]), d; method=method)
        ll_ρ_μ_1_hat[i] =  log_likelihood(parameters(θ_0.π_1, θ_0.π_2, μ_1_hat, θ_0.μ_2, θ_0.σ_1, θ_0.σ_2, ρ[i]), d; method=method)
    end

    # plot log-likelihood
    plot(ρ, ll_ρ_μ_1_0, label = "Log-likelihood for μ_1 = μ_1_0")
    plot!(ρ, ll_ρ_μ_1_hat, label = "Log-likelihood for μ_1 = μ_1_hat")
    vline!([ρ_hat], label = "ρ_hat")
    vline!([θ_0.ρ], label = "ρ_0", legend = :bottomleft)
    xlabel!("ρ")
    ylabel!("Log Likelihood")
end

function fraction_choosing_1(θ::parameters)
    cdf(Normal(0, sqrt(θ.σ_1^2 + θ.σ_2^2 - 2 * θ.ρ * θ.σ_1 * θ.σ_2)), log(θ.π_1) + θ.μ_1 - log(θ.π_2) - θ.μ_2 )
end

# create model_fit table
function model_fit(μ_1_hat::Float64, ρ_hat::Float64, μ_1_hat_s::Float64, ρ_hat_s::Float64, θ_0::parameters, d::data)

    # description of the moments we're looking
    description = ["Fraction choosing occupation 1", "Average observed W1", "Average observed W2", "Std of observed W1", "Std of observed W2"]
    
    # using the true parameters
    s_0 = simulate_data(θ_0; seed =123, N = 100000000)
    W_1 = s_0.W[Bool.(s_0.D)]
    W_2 = s_0.W[Bool.(1 .-s_0.D)]
    true_moments = round.([fraction_choosing_1(θ_0), mean(W_1), mean(W_2), std(W_1), std(W_2)], digits = 3)

    # in the data
    W_1 = d.W[Bool.(d.D)]
    W_2 = d.W[Bool.(1 .-d.D)]

    data_moments = round.([mean(d.D), mean(W_1), mean(W_2), std(W_1), std(W_2)], digits = 3)

    # using mle estimates
    θ_hat = parameters(θ_0.π_1, θ_0.π_2, μ_1_hat, θ_0.μ_2, θ_0.σ_1, θ_0.σ_2, ρ_hat)
    s_hat = simulate_data(θ_hat; seed = 123, N = 100000000)
    W_1 = s_hat.W[Bool.(s_hat.D)]
    W_2 = s_hat.W[Bool.(1 .-s_hat.D)]
    mle_moments = round.([fraction_choosing_1(θ_hat),mean(W_1), mean(W_2), std(W_1), std(W_2)], digits = 3)

    # using the simulation-based mle estimates
    θ_hat_s = parameters(θ_0.π_1, θ_0.π_2, μ_1_hat_s, θ_0.μ_2, θ_0.σ_1, θ_0.σ_2, ρ_hat_s)
    s_hat_s = simulate_data(θ_hat_s; seed = 123, N = 100000000)
    W_1 = s_hat_s.W[Bool.(s_hat_s.D)]
    W_2 = s_hat_s.W[Bool.(1 .-s_hat_s.D)]
    mle_s_moments = round.([fraction_choosing_1(θ_hat_s),mean(W_1), mean(W_2), std(W_1), std(W_2)], digits = 3)
    
    DataFrame(description = description , truth = true_moments, data = data_moments, mle = mle_moments, simulation =mle_s_moments)

end


# create model_fit table
function minimum_wage_counterfactual(μ_1_hat::Float64, ρ_hat::Float64, μ_1_hat_s::Float64, ρ_hat_s::Float64, θ_0::parameters)

    θ_hat = parameters(θ_0.π_1, θ_0.π_2, μ_1_hat, θ_0.μ_2, θ_0.σ_1, θ_0.σ_2, ρ_hat)
    s = simulate_data(θ_hat; seed = 123, N = 10000000)

    W_1_bar = 0.0:0.1:12.0
    frac_choosing_1 = zeros(length(W_1_bar))
    average_W_1 = zeros(length(W_1_bar))
    average_W_2 = zeros(length(W_1_bar))
    std_W_1 = zeros(length(W_1_bar))
    std_W_2 = zeros(length(W_1_bar))

    for i = 1:length(W_1_bar)

        W_1_tilde = max.(s.W_1, W_1_bar[i])

        D_tilde = W_1_tilde .>= s.W_2

        W_tilde = D_tilde .*  W_1_tilde  .+ (1 .- D_tilde) .* s.W_2

        W_1 = W_tilde[Bool.(D_tilde)]
        W_2 = W_tilde[Bool.(1 .- D_tilde)]

        frac_choosing_1[i] = mean(D_tilde)
        average_W_1[i] = mean(W_1)
        average_W_2[i] = mean(W_2)
        std_W_1[i] = std(W_1)
        std_W_2[i] = std(W_2)
    end

    DataFrame(W_1_bar = W_1_bar , frac_choosing_1 = frac_choosing_1, average_W_1 = average_W_1, average_W_2 = average_W_2, std_W_1 =std_W_1,  std_W_2 =std_W_2)

end