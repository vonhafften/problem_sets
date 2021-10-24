##########################################################################################
# Problem set 7
# ECON 899A: Computational Economics
#
# Alex von Hafften
# October 27, 2021
##########################################################################################

using Plots, Tables, DataFrames, CSV

include("model.jl")

D = Initialize_True_Data(; seed = 1234);

# Plot true data
plot(D.x_0, legend = false)
savefig("figures/q2.png");

# smm using different moments
results_4 = estimate(1:2; seed_0 = 1234, seed_1 = 2345)
results_5 = estimate(2:3; seed_0 = 1234, seed_1 = 2345)
results_6 = estimate(1:3; seed_0 = 1234, seed_1 = 2345)

# objective function surface for first stage
plot_J_TH_surface(results_4, D, 1)
savefig("figures/q4a.png");
plot_J_TH_surface(results_5, D, 1)
savefig("figures/q5a.png");
plot_J_TH_surface(results_6, D, 1)
savefig("figures/q6a.png");

# objective function surface for second stage
plot_J_TH_surface(results_4, D, 2)
savefig("figures/q4b.png");
plot_J_TH_surface(results_5, D, 2)
savefig("figures/q5b.png");
plot_J_TH_surface(results_6, D, 2)
savefig("figures/q6b.png");

# Save jacobians
CSV.write("tables/jacobian_4.csv", Tables.table(results_4.jacobian))
CSV.write("tables/jacobian_5.csv", Tables.table(results_5.jacobian))
CSV.write("tables/jacobian_6.csv", Tables.table(results_6.jacobian))

# Summary table
summary_table = create_table([results_4, results_5, results_6], D)
CSV.write("tables/summary.csv", summary_table)

##########################################################################################
# Question 6e
##########################################################################################

results_6_bootstrap = bootstrap_se(1:3)

histogram(results_6_bootstrap[:, 1], fillalpha=0.5, label = "ρ_hat_1");
histogram!(results_6_bootstrap[:, 2], fillalpha=0.5, label = "ρ_hat_2")

savefig("figures/q6e_rho")

histogram(results_6_bootstrap[:, 3], fillalpha=0.5, label = "σ_hat_1");
histogram!(results_6_bootstrap[:, 4], fillalpha=0.5, label = "σ_hat_2")

savefig("figures/q6e_sigma")