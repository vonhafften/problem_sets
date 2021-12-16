####################################################################################################################

# FIN 971 Corporate Finance
# Problem Set 5

# Gomes (2001)
# Alex von Hafften
# December 15, 2021

# This file runs the model and creates charts.

####################################################################################################################

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/fin_971/ps5");

include("model.jl")

@time R_w_friction = Solve_model(0.08, 0.28)
@time R_wo_friction = Solve_model(0.0, 0.0)

####################################################################################################################
# Charts
####################################################################################################################

using Plots
G = Grids()

####################################################################################################################

# Franchise value plot
plot(G.k_grid, R_w_friction.vf, label = false, color="lightgray");
plot!(G.k_grid, R_w_friction.vf[:,[1,5,10]], label = ["z_L" "z_M" "z_H"], legend = :topleft, color = ["green" "red" "blue"]);

savefig("value_function.png")

####################################################################################################################

# Capital policy function
p1 = plot(G.k_grid, R_w_friction.k_pf, label = false, color="lightgray");
plot!(G.k_grid, R_w_friction.k_pf[:,[1,5,10]], label = ["z_L" "z_M" "z_H"], color = ["green" "red" "blue"]);
plot!(G.k_grid, G.k_grid, label = "45 degree", legend = false, color = "black", title = "Capital");

# Investment policy function
p2 = plot(G.k_grid, R_w_friction.i_pf, label = false, color="lightgray", title = "Investment");
plot!(G.k_grid, R_w_friction.i_pf[:,[1,5,10]], label = ["z_L" "z_M" "z_H"], legend = false, color = ["green" "red" "blue"]);

# Exit policy function
p3 = plot(G.k_grid, R_w_friction.x_pf, label = false, color="lightgray", title = "Continuation");
plot!(G.k_grid, R_w_friction.x_pf[:,[1,5,10]], label = ["z_L" "z_M" "z_H"], legend = :bottomright, color = ["green" "red" "blue"]);

# Dividend policy function
p4 = plot(G.k_grid, R_w_friction.d_pf, label = false, color="lightgray", title = "Dividend");
plot!(G.k_grid, R_w_friction.d_pf[:,[1,5,10]], label = ["z_L" "z_M" "z_H"], legend = false, color = ["green" "red" "blue"]);

plot(p1, p2, p3, p4, layout =  (2, 2));
savefig("polcy_functions.png")

####################################################################################################################

# CDF
plot(G.k_grid, R_w_friction.cdf, label = false, color="lightgray");
plot!(G.k_grid, R_w_friction.cdf[:,[1,5,10]], label = ["z_L" "z_M" "z_H"], legend = :bottomright, color = ["green" "red" "blue"]);
savefig("cdf.png")
