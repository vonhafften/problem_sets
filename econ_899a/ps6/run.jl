####################################################################################################################
# ECON 899A Computational Economics
# Problem Set 6

# Hopenhayn and Rogerson (1993)
# Alex von Hafften
# October 20, 2021

# This file runs the model.
####################################################################################################################

using Plots, DataFrames

include("model.jl");



standard = Solve_model()
tv1_1 = Solve_model(;α = 1.0)
tv1_2 = Solve_model(;α = 2.0)

plot(standard.x, label = "Standard")
plot!(tv1_1.x, label = "TV1 Shocks α = 1")
plot!(tv1_2.x, label = "TV1 Shocks α = 2")

standard = Solve_model(;c_f = 15.0)
tv1_1 = Solve_model(;c_f = 15.0, α = 1.0)
tv1_2 = Solve_model(;c_f = 15.0, α = 2.0)

plot(standard.x, label = "Standard")
plot!(tv1_1.x, label = "TV1 Shocks α = 1")
plot!(tv1_2.x, label = "TV1 Shocks α = 2")

create_table([standard_c_f_10, tv1_α_1_c_f_10, standard_c_f_15])