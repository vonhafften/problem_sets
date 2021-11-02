##########################################################################################
# Problem set 7
# ECON 899A: Computational Economics
#
# Alex von Hafften
# October 27, 2021
#
# Estimates coefficent of AR(1) using GMM
##########################################################################################

using Plots, Tables, DataFrames, CSV

include("model.jl");

D_seed = 1234;
R_seed = 202;

# smm using different moments
results_4 = estimate(1:2; D_seed = D_seed, R_seed = R_seed)
results_5 = estimate(2:3; D_seed = D_seed, R_seed = R_seed)
results_6 = estimate(1:3; D_seed = D_seed, R_seed = R_seed)

true_data = Initialize_True_Data(; seed = D_seed);

##########################################################################################
# Plot figures
##########################################################################################

# Plot true data
plot(true_data.x_0, legend = false)
savefig("figures/q2.png");

# objective function surface for first stage
plot_J_TH_surface(results_4, true_data, 1)
savefig("figures/q4a.png");
plot_J_TH_surface(results_5, true_data, 1)
savefig("figures/q5a.png");
plot_J_TH_surface(results_6, true_data, 1)
savefig("figures/q6a.png");

# objective function surface for second stage
plot_J_TH_surface(results_4, true_data, 2)
savefig("figures/q4b.png");
plot_J_TH_surface(results_5, true_data, 2)
savefig("figures/q5b.png");
plot_J_TH_surface(results_6, true_data, 2)
savefig("figures/q6b.png");

##########################################################################################
# Save tables
##########################################################################################

# Save jacobians
CSV.write("tables/jacobian_1_4.csv", Tables.table(results_4.jacobian_1));
CSV.write("tables/jacobian_2_4.csv", Tables.table(results_4.jacobian_2));
CSV.write("tables/jacobian_1_5.csv", Tables.table(results_5.jacobian_1));
CSV.write("tables/jacobian_2_5.csv", Tables.table(results_5.jacobian_2));
CSV.write("tables/jacobian_1_6.csv", Tables.table(results_6.jacobian_1));
CSV.write("tables/jacobian_2_6.csv", Tables.table(results_6.jacobian_2));

# Summary table
summary_table = create_table([results_4, results_5, results_6])
CSV.write("tables/summary.csv", summary_table);

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