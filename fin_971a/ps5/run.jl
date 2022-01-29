####################################################################################################################

# FIN 971 Corporate Finance
# Problem Set 5

# Gomes (2001)
# Alex von Hafften
# December 15, 2021

# This file runs the model and creates charts.

####################################################################################################################

using Plots, DataFrames, Tables, CSV, RegressionTables

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/fin_971/ps5");

include("model.jl")

G = Grids()

@time R = Solve_model(0.08, 0.28)
@time R_frictionless = Solve_model(0.0, 0.0)

####################################################################################################################
# Charts
####################################################################################################################

table_1 = DataFrame(Variable = ["Wage", "Entrant Mass", "Output", "Labor", "Profits"],
                    Baseline = round.([R.wage, R.B, R.Y, R.L_d, R.Π], digits = 4),
                    Frictionless = round.([R_frictionless.wage, R_frictionless.B, R_frictionless.Y, R_frictionless.L_d, R_frictionless.Π], digits = 4))

table_2 = DataFrame(Variable = ["Investment Share", "Share of Financing Costs", "External Financing Cost/Total Costs", "Floatation costs/External Financing Costs"], 
                    Model = round.([R.I/R.Y, R.Λ/R.Y, R.Λ/(R.Λ + R.production_cost), R.floatation_cost/R.Λ], digits = 4))

table_3 = DataFrame(Variable = ["Mean size (capital)", "Mean Investment Rate", "Std. Dev. Investment Rate", "Mean Tobin Q", "Average Cash Flow over Capital", "Std. Dev. Cash Flow over Capital", "Frac. Negative Investment"],
                    Model = round.([R.size_mean, R.i_k_mean, R.i_k_std, R.tobin_q, R.cash_flow_mean, R.cash_flow_std, R.neg_i_frac], digits = 4))

table_4 = DataFrame(Type = ["Receiving seasoned equity \$(d < 0)\$", "Constrained \$(d = 0)\$", "Unconstrained \$(d > 0)\$"],
                    Fraction = round.(R.breakdown_constrained, digits = 4))

CSV.write("table_1.csv", table_1)
CSV.write("table_2.csv", table_2)
CSV.write("table_3.csv", table_3)
CSV.write("table_4.csv", table_4)

####################################################################################################################
i_z_to_plot = [1,6,10]

# Franchise value plot
plot(G.k_grid, R.vf, label = false, color="lightgray");
plot!(G.k_grid, R.vf[:,i_z_to_plot], label = ["z_L" "z_M" "z_H"], legend = :topleft, color = ["green" "red" "blue"])

savefig("value_function.png")

plot(G.k_grid, R_frictionless.vf, label = false, color="lightgray");
plot!(G.k_grid, R_frictionless.vf[:,i_z_to_plot], label = ["z_L" "z_M" "z_H"], legend = :topleft, color = ["green" "red" "blue"])

savefig("value_function_frictionless.png")

####################################################################################################################

# Capital policy function
p1 = plot(G.k_grid, R.k_pf, label = false, color="lightgray");
plot!(G.k_grid, R.k_pf[:,i_z_to_plot], label = ["z_L" "z_M" "z_H"], color = ["green" "red" "blue"]);
plot!(G.k_grid, G.k_grid, label = "45 degree", legend = false, color = "black", title = "Capital");

# Investment policy function
p2 = plot(G.k_grid, R.i_pf, label = false, color="lightgray", title = "Investment");
plot!(G.k_grid, R.i_pf[:,i_z_to_plot], label = ["z_L" "z_M" "z_H"], legend = false, color = ["green" "red" "blue"]);

# Exit policy function
p3 = plot(G.k_grid, R.x_pf, label = false, color="lightgray", title = "Continuation");
plot!(G.k_grid, R.x_pf[:,i_z_to_plot], label = ["z_L" "z_M" "z_H"], legend = :bottomright, color = ["green" "red" "blue"]);

# Dividend policy function
p4 = plot(G.k_grid, R.d_pf, label = false, color="lightgray", title = "Dividend");
plot!(G.k_grid, R.d_pf[:,i_z_to_plot], label = ["z_L" "z_M" "z_H"], legend = false, color = ["green" "red" "blue"]);

plot(p1, p2, p3, p4, layout =  (2, 2))
savefig("policy_functions.png")

####################################################################################################################

# Capital policy function
p1 = plot(G.k_grid, R_frictionless.k_pf, label = false, color="lightgray");
plot!(G.k_grid, R_frictionless.k_pf[:,i_z_to_plot], label = ["z_L" "z_M" "z_H"], color = ["green" "red" "blue"]);
plot!(G.k_grid, G.k_grid, label = "45 degree", legend = false, color = "black", title = "Capital");

# Investment policy function
p2 = plot(G.k_grid, R_frictionless.i_pf, label = false, color="lightgray", title = "Investment");
plot!(G.k_grid, R_frictionless.i_pf[:,i_z_to_plot], label = ["z_L" "z_M" "z_H"], legend = false, color = ["green" "red" "blue"]);

# Exit policy function
p3 = plot(G.k_grid, R_frictionless.x_pf, label = false, color="lightgray", title = "Continuation");
plot!(G.k_grid, R_frictionless.x_pf[:,i_z_to_plot], label = ["z_L" "z_M" "z_H"], legend = :bottomright, color = ["green" "red" "blue"]);

# Dividend policy function
p4 = plot(G.k_grid, R_frictionless.d_pf, label = false, color="lightgray", title = "Dividend");
plot!(G.k_grid, R_frictionless.d_pf[:,i_z_to_plot], label = ["z_L" "z_M" "z_H"], legend = false, color = ["green" "red" "blue"]);

plot(p1, p2, p3, p4, layout =  (2, 2))
savefig("policy_functions_frictionless.png")


####################################################################################################################

# CDF
plot(G.k_grid, R.cdf, label = false, color="lightgray");
plot!(G.k_grid, R.cdf[:,i_z_to_plot], label = ["z_L" "z_M" "z_H"], legend = :bottomright, color = ["green" "red" "blue"]);
plot!(G.k_grid,cumsum(sum(R.μ; dims = 2), dims = 1)/ sum(R.μ), color = "black", label = "Unconditional")
savefig("cdf.png")


plot(G.k_grid, R_frictionless.cdf, label = false, color="lightgray");
plot!(G.k_grid, R_frictionless.cdf[:,i_z_to_plot], label = ["z_L" "z_M" "z_H"], legend = :bottomright, color = ["green" "red" "blue"]);
plot!(G.k_grid,cumsum(sum(R_frictionless.μ; dims = 2), dims = 1)/ sum(R_frictionless.μ), color = "black", label = "Unconditional")
savefig("cdf_frictionless.png")

####################################################################################################################

states = simulate_shocks()

S = simulate_model(R, states)
fhp_baseline = estimate_fhp_regression(S)

S_wo_friction = simulate_model(R_wo_friction, states)
fhp_wo_friction = estimate_fhp_regression(S_wo_friction)

regtable(fhp_baseline, fhp_wo_friction; renderSettings = latexOutput("fhp_table.tex"))