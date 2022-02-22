# Alex von Hafften 
# FIN 970: Asset Pricing
# HW 1 Problem 2
# Professor Ivan Shaliastovich

# This code runs Gibbs sampling to estimate AR(1) process from ./gibbs.jl

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/fin_970/ps1")

include("gibbs.jl")

using Plots, StatsBase

M_1 = Initialize()
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

# scatterplots with lags

scatter(M_1.mu[1:end-1], M_1.rho[2:end], legend = false, title = "Scatterplot of μ and ρ lag", smooth=true)
xlabel!("μ")
ylabel!("ρ_lag")
savefig("p2_q3_mu_rho_lag.png")


# question 4 - posterior mean and standard deviations

mean(M_1.mu)
std(M_1.mu)

mean(M_1.rho)
std(M_1.rho)

mean(M_1.sigma)
std(M_1.sigma)


# question 5 - truncated rho distributions

# use truncate_normal distribution
M_2 = Initialize()
Simulate_MCMC!(M_2; rho_distribution = "truncated_normal")

plot(M_2.rho, legend = false, title = "MCMC of ρ")
ylabel!("ρ")
savefig("p2_q5_rho.png")

# use independence metropolis-hastings
M_3 = Initialize()
Simulate_MCMC!(M_3; rho_distribution = "imh")

plot(M_3.rho, legend = false, title = "MCMC of ρ")
ylabel!("ρ")
savefig("p2_q5_rho.png")
