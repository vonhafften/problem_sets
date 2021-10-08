# Conesa and Krueger (1999)
# Alex von Hafften
# October 6, 2021

# ECON 899A Computational Economics
# Problem Set 4

# This file calls transition.jl and runs the computations

using Plots

include("transition.jl");

# Initial and terminal steady states
θ_0 = 0.11
θ_1 = 0.0

# Guesses for initial steady state
k_guess_0 = 3.513795391470152
l_guess_0 = 0.3468221868691036

# Guesses for terminal steady state
l_guess_1 = 0.36525017624196743
k_guess_1 = 4.598758704122076

################################################################################
######################## Exercise 1 ############################################
################################################################################

@elapsed exercise_1 = Solve_transition(θ_0, θ_1, k_guess_0, k_guess_1, l_guess_0, l_guess_1;
                                       progress = true)

############################ Plots #############################################

plot([exercise_1.r_path repeat([exercise_1.ss_0.r], exercise_1.N_t) repeat([exercise_1.ss_1.r], exercise_1.N_t)],
     label = ["Transition" "Initial SS" "Terminal SS"],
     title = "Interest Rates")

savefig("figures/exercise_1_r.png")

plot([exercise_1.w_path repeat([exercise_1.ss_0.w], exercise_1.N_t) repeat([exercise_1.ss_1.w], exercise_1.N_t)],
     label = ["Transition" "Initial SS" "Terminal SS"],
     title = "Wages",
     legend = :bottomright)

savefig("figures/exercise_1_w.png")

plot([exercise_1.k_demand_path repeat([exercise_1.ss_0.k], exercise_1.N_t) repeat([exercise_1.ss_1.k], exercise_1.N_t)],
     label = ["Transition" "Initial SS" "Terminal SS"],
     title = "Capital",
     legend = :bottomright)

savefig("figures/exercise_1_k.png")

plot([exercise_1.l_demand_path repeat([exercise_1.ss_0.l], exercise_1.N_t) repeat([exercise_1.ss_1.l], exercise_1.N_t)],
     label = ["Transition" "Initial SS" "Terminal SS"],
     title = "Effective Labor",
     legend = :bottomright)

savefig("figures/exercise_1_l.png")

plot(reshape(sum(exercise_1.ev .* exercise_1.ss_0.μ, dims = (2, 3)), 66),
     title = "Average Consumption EV by Age",
     legend = false)

savefig("figures/exercise_1_ev.png")

exercise_1.vote_share

################################################################################
######################## Exercise 2 ############################################
################################################################################

@elapsed exercise_2 = Solve_transition(θ_0, θ_1, k_guess_0, k_guess_1, l_guess_0, l_guess_1;
                                       progress = true, implementation_date = 21, N_t = 50)

############################ Plots #############################################

plot(exercise_2.θ,
     label = "Transition",
     title = "Theta")

savefig("figures/exercise_2_theta.png")

plot([exercise_2.b_path repeat([exercise_2.ss_0.b], exercise_2.N_t) repeat([exercise_2.ss_1.b], exercise_2.N_t)],
     label = ["Transition" "Initial SS" "Terminal SS"],
     title = "Social Security Benefit")

savefig("figures/exercise_2_b.png")

plot([exercise_2.r_path repeat([exercise_2.ss_0.r], exercise_2.N_t) repeat([exercise_2.ss_1.r], exercise_2.N_t)],
     label = ["Transition" "Initial SS" "Terminal SS"],
     title = "Interest Rates")

savefig("figures/exercise_2_r.png")

plot([exercise_2.w_path repeat([exercise_2.ss_0.w], exercise_2.N_t) repeat([exercise_2.ss_1.w], exercise_2.N_t)],
     label = ["Transition" "Initial SS" "Terminal SS"],
     title = "Wages",
     legend = :bottomright)

savefig("figures/exercise_2_w.png")

plot([exercise_2.k_demand_path repeat([exercise_2.ss_0.k], exercise_2.N_t) repeat([exercise_2.ss_1.k], exercise_2.N_t)],
     label = ["Transition" "Initial SS" "Terminal SS"],
     title = "Capital",
     legend = :bottomright)

savefig("figures/exercise_2_k.png")

plot([exercise_2.l_demand_path repeat([exercise_2.ss_0.l], exercise_2.N_t) repeat([exercise_2.ss_1.l], exercise_2.N_t)],
     label = ["Transition" "Initial SS" "Terminal SS"],
     title = "Effective Labor",
     legend = :bottomright)

savefig("figures/exercise_2_l.png")

plot([reshape(sum(exercise_2.ev .* exercise_2.ss_0.μ, dims = (2, 3)), 66) reshape(sum(exercise_1.ev .* exercise_1.ss_0.μ, dims = (2, 3)), 66)],
     title = "Average Consumption EV by Age",
     label = ["Expected" "Unexpected"])

savefig("figures/exercise_2_ev.png")

exercise_2.vote_share
