# Alex von Hafften 
# FIN 970: Asset Pricing
# HW 1 Problem 2
# Professor Ivan Shaliastovich

# This code use Gibbs sampling to estimate AR(1) process.
# Normal drift, normal persistance, inverse-gamma variance.

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/fin_970/ps1")

using XLSX, DataFrames, Parameters, Distributions

@with_kw struct Priors 
    # for drift
    m = 0.3
    s = 0.5

    # for persistance
    rho_tilde = 0.95
    omega     = 0.2

    # for variance
    alpha = 6.0
    beta  = 4.0
end

@with_kw struct Data_Moments
    # data
    data = DataFrame(XLSX.readtable("Long_yield_data.xlsx", "data")...)
    yields = Float64.(data.IGUSA10D)
    y = yields[2:end]
    y_lag = yields[1:(end-1)]
    T = size(y)[1]

    # data moments
    y_2_bar = 1/T * sum(y.^2)
    y_lag_2_bar = 1/T * sum(y_lag.^2)
    y_bar = 1/T * sum(y)
    y_lag_bar = 1/T * sum(y_lag)
    z_bar = 1/T * sum(y .* y_lag)
end 

mutable struct MCMC
    N::Int64                   # length of mcmc
    burn_in::Int64             # number of burn-in periods
    mu::Vector{Float64}        # vector of draws for m
    rho::Vector{Float64}       # vector of draws for rho
    sigma::Vector{Float64}     # vector of draws for omega
end


function Initialize_MCMC()
    P = Priors()

    N = 1000
    burn_in = 100

    mu    = zeros(N+burn_in)
    rho   = zeros(N+burn_in)
    sigma = zeros(N+burn_in)

    mu[1]    = P.m
    rho[1]   = P.rho_tilde
    sigma[1] = P.beta/(P.alpha - 1)

    return MCMC(N, burn_in, mu, rho, sigma)
end

function Simulate_MCMC!(M::MCMC; rho_distribution::String = "normal")
    P = Priors()
    D = Data_Moments()

    for i = 2:(M.N+M.burn_in)

        # draw mu
        nu_mu   = (M.sigma[i-1]^2)/(P.s^2)
        m_tilde = P.m * (nu_mu/(nu_mu + D.T)) + (M.rho[i-1] * D.y_lag_bar-D.y_bar)*(D.T/(nu_mu + D.T))
        s_tilde = sqrt(M.sigma[i-1]^2/(nu_mu + D.T))
        M.mu[i] = rand(Normal(m_tilde, s_tilde))

        # draw rho
        nu_rho          = (M.sigma[i-1]^2)/(P.omega^2)
        rho_tilde_tilde = (nu_rho/(nu_rho + D.T*D.y_lag_2_bar)) * P.rho_tilde +
                        (D.T/(nu_rho + D.T * D.y_lag_2_bar)*(D.z_bar - M.mu[i]*D.y_lag_bar))
        omega_tilde     = sqrt((M.sigma[i-1]^2)/(nu_rho + D.T*D.y_lag_2_bar))
        
        if rho_distribution == "normal"
            M.rho[i] = rand(Normal(rho_tilde_tilde, omega_tilde))
        elseif rho_distribution == "truncated_normal"
            M.rho[i] = rand(truncated(Normal(rho_tilde_tilde, omega_tilde), -1.0, 1.0))
        elseif rho_distribution == "imh"
            proposal = Normal(rho_tilde_tilde, omega_tilde)

            # get new candidate
            rho_candidate = rand(proposal)

            # get f of rho_candidate
            if rho_candidate > 1.0
                f_rho_candidate = 0.0
            elseif rho_candidate < -1.0
                f_rho_candidate = 0.0
            else 
                f_rho_candidate = pdf(proposal, rho_candidate)
            end

            # get f of rho_tilde_tilde
            if rho_tilde_tilde > 1.0
                f_rho_tilde_tilde = 0.0
            elseif rho_tilde_tilde < -1.0
                f_rho_tilde_tilde = 0.0
            else 
                f_rho_tilde_tilde = pdf(proposal, rho_tilde_tilde)
            end

            q_rho_candidate = pdf(proposal, rho_candidate)
            q_rho_tilde_tilde = pdf(Normal(rho_candidate, omega_tilde), rho_tilde_tilde)
            
            threshold = min((f_rho_candidate/f_rho_tilde_tilde)*(q_rho_tilde_tilde/q_rho_candidate), 1)

            if isnan(threshold)
                if f_rho_candidate == 0.0
                    threshold = 0.0
                else
                    threshold = 1.0
                end
            end

            U = rand()

            if U < threshold
                M.rho[i] = rho_candidate
            else
                M.rho[i] = M.rho[i-1]
            end
        else
            error("specify valid method for rho_distribution: normal, truncate_normal, imh")
        end

        # draw sigma
        alpha_tilde = P.alpha + D.T
        beta_tilde = P.beta + D.T * (D.y_2_bar + M.mu[i]^2 + M.rho[i]^2*D.y_lag_2_bar - 
                            2 * M.mu[i] *D.y_bar - 2 * M.rho[i] *D.z_bar + 2*M.rho[i]*M.mu[i]*D.y_lag_bar)
        M.sigma[i] = rand(InverseGamma(alpha_tilde/2, beta_tilde/2))
    end

    # drop burn-in
    M.mu = M.mu[(M.burn_in+1):(M.N + M.burn_in)]
    M.rho = M.rho[(M.burn_in+1):(M.N + M.burn_in)]
    M.sigma = M.sigma[(M.burn_in+1):(M.N + M.burn_in)]

    return M
end

mutable struct Forecast
    mean::Vector{Float64}
    ci::Matrix{Float64}
end

function forecast(M::MCMC, N_ahead::Int64)
    D = Data_Moments()

    mean_forecast = zeros(D.T)
    ci_forecast = zeros(D.T, 2)

    for t=1:D.T
        forecasts = M.mu .*((1 .- M.rho.^N_ahead)./(1 .- M.rho)) .+ M.rho .^ N_ahead * D.yields[t]
        mean_forecast[t] = mean(forecasts)
        ci_forecast[t, 1] = quantile(forecasts, 0.025)
        ci_forecast[t, 2] = quantile(forecasts, 0.975)
    end
    return Forecast(mean_forecast, ci_forecast)
end