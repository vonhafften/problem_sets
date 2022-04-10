# Alex von Hafften
# Problem set 1
# ECON 717B: Applied Econometrics
# Professor Matt Wiswall
# April 7, 2022

# This file calls functions from ./model.jl

cd("/Users/vonhafften/Documents/UW Madison/problem_sets/econ_717b/ps1/")
include("model.jl")

using Plots, CSV

################################################################################################
# part a - why we can normalize π_1 and π_2
################################################################################################

# see write-up

################################################################################################
# part b - program to compute the model
################################################################################################

# see run.jl

################################################################################################
# part c - choose vector of parameters such that 60 percent of simulations choose occupation 1
################################################################################################

# true parameters
θ_0 = parameters(1.0, 1.0, 1.15, 1.0, 0.4, 0.5, 0.25)

# simulate true data
data_0 = simulate_data(θ_0;seed = 131)

# fraction that choose occupation is about 60 percent
sum(data_0.D)/data_0.N

################################################################################################
# part d - write mle estimator for θ
################################################################################################

# see write-up

################################################################################################
# part e - write code to estimate θ using mle
################################################################################################

(μ_1_hat, ρ_hat) = mle(θ_0, data_0)

(μ_1_hat_s, ρ_hat_s) = mle(θ_0, data_0; method = "simulation")

################################################################################################
# part f - plot identification figure
################################################################################################

id_plot_3d(μ_1_hat, ρ_hat, θ_0, data_0)
title!("MLE Id. Plot: μ_1_hat and ρ_hat")
savefig("q_f_3d.png")

id_plot_μ_1(μ_1_hat, ρ_hat, θ_0, data_0)
title!("MLE Id. Plot: μ_1_hat")
savefig("q_f_mu_1.png")

id_plot_ρ(μ_1_hat, ρ_hat, θ_0, data_0)
title!("MLE Id. Plot: ρ_hat")
savefig("q_f_rho.png")

id_plot_3d(μ_1_hat_s, ρ_hat_s, θ_0, data_0; method = "simulation")
title!("Simulated MLE Id. Plot: μ_1_hat and ρ_hat")
savefig("q_f_3d_s.png")

id_plot_μ_1(μ_1_hat_s, ρ_hat_s, θ_0, data_0;method = "simulation")
title!("Simulated MLE Id. Plot: μ_1_hat")
savefig("q_f_mu_1_s.png")

id_plot_ρ(μ_1_hat_s, ρ_hat_s, θ_0, data_0;method = "simulation")
title!("Simulated MLE Id. Plot: ρ_hat")
savefig("q_f_rho_s.png")

################################################################################################
# part g - table of parameters and model fit
################################################################################################

model_fit_table = model_fit(μ_1_hat, ρ_hat, μ_1_hat_s, ρ_hat_s, θ_0, data_0)

CSV.write("q_g_table_1.csv", select(model_fit_table, :description, :data, :mle))
CSV.write("q_g_table_2.csv", model_fit_table)

################################################################################################
# part h - minimum wage counterfactual
################################################################################################

mw_cf = minimum_wage_counterfactual(μ_1_hat, ρ_hat, μ_1_hat_s, ρ_hat_s, θ_0)

plot(mw_cf.W_1_bar, mw_cf.frac_choosing_1, legend = false)
xlabel!("Minimum Wage in Occupation 1")
ylabel!("Fraction in Occupation 1")
savefig("q_h_fraction.png")

plot(mw_cf.W_1_bar, mw_cf.average_W_1, legend = :topleft, label = "W_1")
plot!(mw_cf.W_1_bar, mw_cf.average_W_2, label = "W_2")
ylabel!("Average Observed Wages")
xlabel!("Minimum Wage in Occupation 1")
savefig("q_h_average.png")

plot(mw_cf.W_1_bar, mw_cf.std_W_1, label = "W_1", legend = :bottomleft)
plot!(mw_cf.W_1_bar, mw_cf.std_W_2, label = "W_2")
ylabel!("Std. Dev. of Observed Wages")
xlabel!("Minimum Wage in Occupation 1")
savefig("q_h_std_dev.png")

CSV.write("q_h_table.csv", round.(mw_cf[[1,30],:], digits=3))