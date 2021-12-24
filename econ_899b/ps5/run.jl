# Alex von Hafften 
# Besanko and Doraszelski
# Supplement to Computational Economics
# December 23, 2021

# This file runs the plotting functions
cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_899b/ps5/")

include("model.jl")
include("plotting.jl")

# run quantity and price competition with 2 depreciation levels

@time model_q_δ1  = Solve_model("quantity", 0.1);
@time model_q_δ3  = Solve_model("quantity", 0.3);

@time model_p_δ1  = Solve_model("price", 0.1);
@time model_p_δ3  = Solve_model("price", 0.3);

# Static results plots
plot_static_results(model_q_δ1);
savefig("q_static.png")

plot_static_results(model_p_δ1);
savefig("p_static.png")

# Dynamic results plots
plot_dynamic_results(model_q_δ1);
savefig("q_low_delta_dynamic.png")

plot_dynamic_results(model_q_δ3);
savefig("q_high_delta_dynamic.png")

plot_dynamic_results(model_p_δ1);
savefig("p_low_delta_dynamic.png")

plot_dynamic_results(model_p_δ3);
savefig("p_high_delta_dynamic.png")