# Problem Set 3 - BLP
# Computational Economics
# Taught by JF Houde

# Alex von Hafften

using DataFrames, StatFiles, LinearAlgebra, Statistics, Plots, Optim

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_899b/ps3")

# read in functions
include("data.jl")
include("model.jl")

####################################################################################
# Question 1 - Invert demand for 1985 and plot sup norm error
####################################################################################

err_list_1985 = invert_demand(1985, 0.6)[2]
err_list_1985_cm = invert_demand(1985, 0.6; tolerence = [1e-12,1e-12])[2]

# plot omitting first iteration with huge error
plot(err_list_1985[2:end], label = "CM + NM");
hline!([1], label="CM + NM Threshold");
plot!(err_list_1985_cm[2:end], label = "CM");
xlabel!("Error");
ylabel!("Iteration #")

savefig("question_1.png")


####################################################################################
# Question 2 - GMM first stage
####################################################################################

λ_grid = 0.0:0.1:1.0

W_1 = inv(Z' * Z) # first stage weighting matrix

# objective function over λ_grid
gmm_obj_1_grid = zeros(length(λ_grid))
for i = 1:length(λ_grid)
    println(i)
    gmm_obj_1_grid[i] = gmm_obj(λ_grid[i], W_1)
end

@time gmm_1 = optimize(λ -> gmm_obj(λ, W_1), 0.0, 1.0, Brent())

plot(λ_grid, gmm_obj_1_grid, label = "First stage objective");
scatter!([gmm_1.minimizer], [gmm_1.minimum], label = "First stage estimate");
xlabel!("λ_p")
savefig("question_2.png")


####################################################################################
# Question 3 - GMM second stage
####################################################################################

ξ_hat = compute_ρ(gmm_1.minimizer, W_1)

# second stage weighting matrix
W_2 = inv((Z .* ξ_hat)' *(Z .* ξ_hat))

# objective function over λ_grid
gmm_obj_2_grid = zeros(length(λ_grid))
for i = 1:length(λ_grid)
    println(i)
    gmm_obj_2_grid[i] = gmm_obj(λ_grid[i], W_2)
end

@time gmm_2 = optimize(λ -> gmm_obj(λ, W_2), 0.0, 1.0, Brent())

plot(λ_grid, gmm_obj_1_grid, label = "1st stage objective");
scatter!([gmm_1.minimizer], [gmm_1.minimum], label = "1st stage estimate");
plot!(λ_grid, gmm_obj_2_grid, label = "2nd stage objective");
scatter!([gmm_2.minimizer], [gmm_2.minimum], label = "2nd stage estimate");
xlabel!("λ_p")

savefig("question_3.png")