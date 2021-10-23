##########################################################################################
# Problem set 7
# ECON 899A: Computational Economics
#
# Alex von Hafften
# October 27, 2021
##########################################################################################

include("model.jl")

##########################################################################################
# Question 2
##########################################################################################

P = Primitives()
plot(P.x_0, legend = false)
savefig("figures/q2.png");

##########################################################################################
# Question 4
##########################################################################################

results_4 = estimate(1:2; seed = 2345)

plot_J_TH_surface(results_4, 1)
savefig("figures/q4a.png");

results_4.ρ_hat_1
results_4.σ_hat_1
results_4.W_hat
results_4.ρ_hat_2
results_4.σ_hat_2

plot_J_TH_surface(results_4, 2)
savefig("figures/q4a.png");

##########################################################################################
# Question 5
##########################################################################################

results_5 = estimate(2:3; seed = 2345)

plot_J_TH_surface(results_5, 1)
savefig("figures/q5a.png");

results_5.ρ_hat_1
results_5.σ_hat_1
results_5.W_hat
results_5.ρ_hat_2
results_5.σ_hat_2

plot_J_TH_surface(results_5, 2)
savefig("figures/q5b.png");

##########################################################################################
# Question 5
##########################################################################################

results_6 = estimate(1:3; seed = 2345)

plot_J_TH_surface(results_6, 1)
savefig("figures/q5a.png");

results_6.ρ_hat_1
results_6.σ_hat_1
results_6.W_hat
results_6.ρ_hat_2
results_6.σ_hat_2

plot_J_TH_surface(results_6, 2)
savefig("figures/q5b.png");