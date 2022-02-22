# Alex von Hafften 
# FIN 970: Asset Pricing
# HW 1 Problem 2
# Professor Ivan Shaliastovich

# This code runs Gibbs sampling to estimate AR(1) process from ./gibbs.jl

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/fin_970/ps1")

include("gibbs.jl")

using Plots, StatsBase

M_1 = Initialize_MCMC()
Simulate_MCMC!(M_1)

# question 3

# line plots

plot(M_1.mu, legend = false, title = "MCMC of μ")
ylabel!("μ")
savefig("p2_q3_mu.png")

plot(M_1.rho, legend = false, title = "MCMC of ρ")
ylabel!("ρ")
savefig("p2_q3_rho.png")

plot(M_1.sigma, legend = false, title = "MCMC of σ")
ylabel!("σ")
savefig("p2_q3_sigma.png")

# ar(1) persistance

# mu
y = M_1.mu[2:end]
x = [ones(M_1.N-1) M_1.mu[1:end-1]]
inv(x'*x)*x'*y

# rho
y = M_1.rho[2:end]
x = [ones(M_1.N-1) M_1.rho[1:end-1]]
inv(x'*x)*x'*y

# sigma
y = M_1.sigma[2:end]
x = [ones(M_1.N-1) M_1.sigma[1:end-1]]
inv(x'*x)*x'*y

# autocorrelation plots

lags = 0:10

acf_mu = autocor(M_1.mu, lags)
plot(lags, acf_mu, line=:stems, marker=:circle, legend = false)
xlabel!("Lag")
ylabel!("Autocorrelation")
title!("ACF μ")
savefig("p2_q3_mu_acf.png")

acf_rho = autocor(M_1.rho, lags)
plot(lags, acf_rho, line=:stems, marker=:circle, legend = false)
xlabel!("Lag")
ylabel!("Autocorrelation")
title!("ACF ρ")
savefig("p2_q3_rho_acf.png")

acf_sigma = autocor(M_1.sigma, lags)
plot(lags, acf_sigma, line=:stems, marker=:circle, legend = false)
xlabel!("Lag")
ylabel!("Autocorrelation")
title!("ACF σ")
savefig("p2_q3_sigma_acf.png")



# scatterplots

scatter(M_1.mu, M_1.rho, legend = false, title = "Scatterplot of μ and ρ", smooth=true)
xlabel!("μ")
ylabel!("ρ")
savefig("p2_q3_mu_rho.png")

scatter(M_1.mu, M_1.sigma, legend = false, title = "Scatterplot of μ and σ", smooth=true)
xlabel!("μ")
ylabel!("σ")
savefig("p2_q3_mu_sigma.png")

scatter(M_1.rho, M_1.sigma, legend = false, title = "Scatterplot of ρ  and σ", smooth=true)
xlabel!("ρ")
ylabel!("σ")
savefig("p2_q3_rho_sigma.png")


# question 4 - posterior mean and standard deviations

mean(M_1.mu)
std(M_1.mu)

mean(M_1.rho)
std(M_1.rho)

mean(M_1.sigma)
std(M_1.sigma)


# question 5/6 - truncated rho distributions

# use truncate_normal distribution
M_2 = Initialize_MCMC()
Simulate_MCMC!(M_2; rho_distribution = "truncated_normal")

plot(M_2.rho, legend = false, title = "MCMC of ρ (Truncated Normal)")
ylabel!("ρ")
savefig("p2_q5_rho_truncated_normal.png")

mean(M_2.mu)
std(M_2.mu)

mean(M_2.rho)
std(M_2.rho)

mean(M_2.sigma)
std(M_2.sigma)

# use independence metropolis-hastings
M_3 = Initialize_MCMC()
Simulate_MCMC!(M_3; rho_distribution = "imh")

plot(M_3.rho, legend = false, title = "MCMC of ρ (IMH)")
ylabel!("ρ")
savefig("p2_q5_rho_imh.png")

mean(M_3.mu)
std(M_3.mu)

mean(M_3.rho)
std(M_3.rho)

mean(M_3.sigma)
std(M_3.sigma)

scatter(M_2.mu, M_2.rho, legend = false, title = "Scatterplot of μ and ρ using truncated normal", smooth=true)
xlabel!("μ")
ylabel!("ρ")
savefig("p2_q6_mu_rho_truncated_normal.png")

scatter(M_3.mu, M_3.rho, legend = false, title = "Scatterplot of μ and ρ using IMH", smooth=true)
xlabel!("μ")
ylabel!("ρ")
savefig("p2_q6_mu_rho_imh.png")

# question 6

forecast_5 = forecast(M_1, 5)

plot(forecast_5.mean, label = "Mean", title = "Five-Year Yield Forecast")
plot!(forecast_5.ci, label = ["CI 2.5%" "CI 97.5%"], line=:dash)
savefig("p2_q6_forecast.png")