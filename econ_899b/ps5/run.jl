# Alex von Hafften 
# Besanko and Doraszelski
# Supplement to Computational Economics
# December 23, 2021

# This file runs the plotting functions

include("model.jl")
include("plotting.jl")

@time model_quantity_competition = Solve_model("quantity", 0.1)
@time model_price_competition = Solve_model("price", 0.1)

plot_static_results(model_quantity_competition)
plot_dynamic_results(model_quantity_competition)

plot_static_results(model_price_competition)
plot_dynamic_results(model_price_competition)
