# Alex von Hafften
# ECON 810: Advanced Macro
# PS 5 - Code for running the model
# Professor Carter Braxton

using Plots

cd("/Users/vonhafften/Documents/uw_madison/problem_sets/econ_810a/ps5/")

include("model.jl")
include("simulation.jl")

P = Primitives()
R = Initialize()

Solve!(R)

# investment at t=8

plot(P.grid_b, R.pf_i_8[:, :, 1], color=:lightgray, label = false, title = "Investment at t = 8 for Child with Lowest HC")
plot!(P.grid_b, R.pf_i_8[:, P.N_h, 1], label = "Parent with Highest HC", color = :red, legend = :topleft)
plot!(P.grid_b, R.pf_i_8[:, 1, 1], label = "Parent with Lowest HC", color = :blue)
xlabel!("Assets")
ylabel!("Investment")

savefig("i_c_l.png")

plot(P.grid_b, R.pf_i_8[:, :, P.N_h], color=:lightgray, label = false, title = "Investment at t = 8 for Child with Highest HC")
plot!(P.grid_b, R.pf_i_8[:, P.N_h, P.N_h], label = "Parent with Highest HC", color = :red, legend = :topleft)
plot!(P.grid_b, R.pf_i_8[:, 1, P.N_h], label = "Parent with Lowest HC", color = :blue)
xlabel!("Assets")
ylabel!("Investment")

savefig("i_c_h.png")

plot(P.grid_b, R.pf_i_8[:, 1, :], color=:lightgray, label = false, title = "Investment at t = 8 by Parent with Lowest HC")
plot!(P.grid_b, R.pf_i_8[:, 1, P.N_h], label = "Child with Highest HC", color = :red, legend = :topleft)
plot!(P.grid_b, R.pf_i_8[:, 1, 1], label = "Child with Lowest HC", color = :blue)
xlabel!("Assets")
ylabel!("Investment")

savefig("i_p_l.png")

plot(P.grid_b, R.pf_i_8[:, P.N_h, :], color=:lightgray, label = false, title = "Investment at t = 8 by Parent with Highest HC")
plot!(P.grid_b, R.pf_i_8[:, P.N_h, P.N_h], label = "Child with Highest HC", color = :red, legend = :topleft)
plot!(P.grid_b, R.pf_i_8[:, P.N_h, 1], label = "Child with Lowest HC", color = :blue)
xlabel!("Assets")
ylabel!("Investment")

savefig("i_p_h.png")

# transfers plots

plot(P.grid_b, R.pf_τ_9[:, :, 1], color=:lightgray, label = false, title = "Transfer for Child with Lowest HC")
plot!(P.grid_b, R.pf_τ_9[:, P.N_h, 1], label = "Parent with Highest HC", color = :red, legend = :bottomright)
plot!(P.grid_b, R.pf_τ_9[:, 1, 1], label = "Parent with Lowest HC", color = :blue)
xlabel!("Assets")
ylabel!("Transfer")

savefig("tau_c_l.png")

plot(P.grid_b, R.pf_τ_9[:, :, P.N_h], color=:lightgray, label = false, title = "Transfer for Child with Highest HC")
plot!(P.grid_b, R.pf_τ_9[:, P.N_h, P.N_h], label = "Parent with Highest HC", color = :red, legend = :bottomright)
plot!(P.grid_b, R.pf_τ_9[:, 1, P.N_h], label = "Parent with Lowest HC", color = :blue)
xlabel!("Assets")
ylabel!("Transfer")

savefig("tau_c_h.png")

plot(P.grid_b, R.pf_τ_9[:, 1, :], color=:lightgray, label = false, title = "Transfer from Parent with Lowest HC")
plot!(P.grid_b, R.pf_τ_9[:, 1, P.N_h], label = "Child with Highest HC", color = :red, legend = :bottomright)
plot!(P.grid_b, R.pf_τ_9[:, 1, 1], label = "Child with Lowest HC", color = :blue, legend = :bottomright)
xlabel!("Assets")
ylabel!("Transfer")

savefig("tau_p_l.png")

plot(P.grid_b, R.pf_τ_9[:, P.N_h, :], color=:lightgray, label = false, title = "Transfer from Parent with Highest HC")
plot!(P.grid_b, R.pf_τ_9[:, P.N_h, P.N_h], label = "Child with Highest HC", color = :red, legend = :bottomright)
plot!(P.grid_b, R.pf_τ_9[:, P.N_h, 1], label = "Child with Lowest HC", color = :blue, legend = :bottomright)
xlabel!("Assets")
ylabel!("Transfer")

savefig("tau_p_h.png")


S = Initialize_simulation(10000)
Simulate!(R, S)

plot(4:12, mean(S.h, dims = 1)', title = "Mean Human Capital", legend = false)
xlabel!("Model Age")
ylabel!("Human Capital")

savefig("mean_hc.png")

plot(4:12, var(S.h, dims = 1)', title = "Cross Sectional Variance of Human Capital", legend = false)
xlabel!("Model Age")
ylabel!("Human Capital^2")

savefig("var_hc.png")

plot(4:12, mean(S.b, dims = 1)', title = "Mean Assets", legend = false)
xlabel!("Model Age")
ylabel!("Assets")

savefig("mean_assets.png")

plot(4:12, var(S.b, dims = 1)', title = "Cross Sectional Variance of Assets", legend = false)
xlabel!("Model Age")
ylabel!("Assets^2")

savefig("var_assets.png")

histogram(S.h_c[:, end], legend = false, title = "Human Capital at t = 4")
xlabel!("Human Capital")
ylabel!("Frequency")

savefig("hc_histogram.png")

histogram(S.τ, legend = false, title = "Transfer")
xlabel!("Transfer")
ylabel!("Frequency")

savefig("transfer_histogram.png")


