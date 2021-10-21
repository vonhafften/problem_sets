####################################################################################################################
# ECON 899A Computational Economics
# Problem Set 6

# Hopenhayn and Rogerson (1993)
# Alex von Hafften
# October 20, 2021

# This file runs the model.
####################################################################################################################

using Plots

include("model.jl");

# Estimate models with c_f = 10.0
standard_10 = Solve_model()
tv1_1_10 = Solve_model(;α = 1.0)
tv1_2_10 = Solve_model(;α = 2.0)

# Estimate models with c_f = 15.0
standard_15 = Solve_model(;c_f = 15.0)
tv1_1_15 = Solve_model(;c_f = 15.0, α = 1.0)
tv1_2_15 = Solve_model(;c_f = 15.0, α = 2.0)

# create table of results
table = create_table([standard_10, tv1_1_10, tv1_2_10, standard_15, tv1_1_15, tv1_2_15])
CSV.write("tables/summary.csv", table)

# Plot results for c_f = 10.0
plot(standard_10.x, label = "Standard")
plot!(tv1_1_10.x, label = "TV1 Shocks α = 1")
plot!(tv1_2_10.x, label = "TV1 Shocks α = 2")
plot!(title = "Exit Decisions for c_f = 10", xlab = "Productivity State #", ylab = "Pr(Exit)")

savefig("figures/c_f_10.png")

# Plot results for c_f = 15.0
plot(standard_15.x, label = "Standard", title = "Exit Decisions")
plot!(tv1_1_15.x, label = "TV1 Shocks α = 1")
plot!(tv1_2_15.x, label = "TV1 Shocks α = 2")
plot!(title = "Exit Decisions for c_f = 10", xlab = "Productivity State #", ylab = "Pr(Exit)")

savefig("figures/c_f_15.png")
