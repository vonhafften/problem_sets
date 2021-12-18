####################################################################################################################

# ECON 899B Computational Economics
# Professor: JF Houde
# Problem Set 4

# Inventory Problem
# Alex von Hafften
# December 17, 2021

# This file contains the model code.

####################################################################################################################

####################################################################################################################
################################## Define primitives and results structures ########################################
####################################################################################################################

using CSV, Tables, DataFrames, LinearAlgebra, Optim, Plots

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/econ_899b/ps4/")

# Model parameters
λ        = -4.0                 # Stockout penalty
α        = 2.0                  # Coefficient on consumption shock
β        = 0.99                 # Discount rate
n        = 36

#########################################################################################
# read in data
state_space_df    = select(DataFrame(CSV.File("PS4/PS4_state_space.csv")), Not(:Column1))        # state space definitions df
transition_a0_df  = select(DataFrame(CSV.File("PS4/PS4_transition_a0.csv")), Not(:Column1))      # transition matrix if a = 0 df
transition_a1_df  = select(DataFrame(CSV.File("PS4/PS4_transition_a1.csv")), Not(:Column1))       # transition matrix if a = 1 df

F_0 = Array(select(transition_a0_df, Not(:id))) # transition matrix if a = 0 matrix
F_1 = Array(select(transition_a1_df, Not(:id))) # transition matrix if a = 1 matrix

#########################################################################################
# Question 1
#########################################################################################

# current payoffs not including T1EV shock
U_1 = α .* state_space_df.C - state_space_df.P
U_0 = (state_space_df.I .> 0) .* α .* state_space_df.C .+ (state_space_df.I .== 0) .* (state_space_df.C .> 0) .* λ

# VFI on expected value
EV = log.(exp.(U_0) .+ exp.(U_1)) .+ Base.MathConstants.eulergamma
i, err, maxiter = 0, 100, 1e+14
while (err > 1e-14) & (i < maxiter)
    i += 1
    EV_next = log.(exp.(U_0 .+ β * F_0 * EV) .+ exp.(U_1 .+ β * F_1 * EV)) .+ Base.MathConstants.eulergamma
    err = maximum(abs.(EV .- EV_next))
    EV = EV_next
    # println("Iteration #", i)
    # println("Current error:", err)
end
print(i)

EV_VFI = EV

table_1 = select(state_space_df, Not(:id))
table_1.U0 = round.(U_0; digits = 3)
table_1.U1 = round.(U_1; digits = 3)
table_1.EV = round.(EV; digits = 3)

CSV.write("table_1.csv", table_1)

#########################################################################################
# Question 2
#########################################################################################

simdata = select(DataFrame(CSV.File("PS4/PS4_simdata.csv")), Not(:Column1))        # state space definitions df

# group by state_id and calculate estimated frequency
estimated_frequencies = combine(groupby(simdata, :state_id), :choice => mean => :P_1)
estimated_frequencies.P_0 = 1 .- estimated_frequencies.P_1

# constrain between 0.0001 and 0.999
estimated_frequencies.P_0 = max.(estimated_frequencies.P_0, 0.001)
estimated_frequencies.P_1 = max.(estimated_frequencies.P_1, 0.001)
estimated_frequencies.P_0 = min.(estimated_frequencies.P_0, 0.999)
estimated_frequencies.P_1 = min.(estimated_frequencies.P_1, 0.999)

# check all between 0.0001 and 0.999
sum(estimated_frequencies.P_0 .> 0.999)
sum(estimated_frequencies.P_1 .> 0.999)
sum(estimated_frequencies.P_0 .< 0.001)
sum(estimated_frequencies.P_1 .< 0.001)

# expectation of T1EV shocks
e_0 = Base.MathConstants.eulergamma .- log.(estimated_frequencies.P_0)
e_1 = Base.MathConstants.eulergamma .- log.(estimated_frequencies.P_1)

# Initial value for CCP
P = (1 .+ exp.(-(estimated_frequencies.P_1 .- estimated_frequencies.P_0))).^(-1)

i, err, maxiter = 0, 100, 1e+14
while (err > 1e-14) & (i < maxiter)
    i += 1
    F = (1 .- P) .* F_0 .+ P .* F_1
    EV = inv(I - β * F) * ((1 .- P).*(U_0 + e_0) + P .* (U_1 + e_1))
    V_tilde = (U_1 .+ β .* F_1 * EV) - (U_0 .+ β .* F_0 * EV)
    P_next = (1 .+ exp.(-V_tilde)).^(-1)
    err = maximum(abs.(P .- P_next))
    P = P_next
    # println("Iteration #", i)
    # println("Current error:", err)
end
print(i)

EV_PFI = EV

table_2 = select(state_space_df, Not(:id))
table_2.P0hat = round.(estimated_frequencies.P_0; digits = 3)
table_2.P1hat = round.(estimated_frequencies.P_1; digits = 3)
table_2.EV = table_1.EV
table_2.EVhat = round.(EV; digits = 3)

CSV.write("table_2.csv", table_2)

maximum(abs.((EV_VFI .- EV_PFI)./EV_VFI))
mean(abs.((EV_VFI .- EV_PFI)./EV_VFI))

#########################################################################################
# Question 4
#########################################################################################

function log_likelihood(λ_temp, simdata, estimated_frequencies, state_space_df, F_0, F_1)
    P = (1 .+ exp.(-(estimated_frequencies.P_1 .- estimated_frequencies.P_0))).^(-1)

    U_1 = α .* state_space_df.C - state_space_df.P
    U_0 = (state_space_df.I .> 0) .* α .* state_space_df.C .+ (state_space_df.I .== 0) .* (state_space_df.C .> 0) .* λ_temp

    i, err, maxiter = 0, 100, 1e+14
    while (err > 1e-14) & (i < maxiter)
        i += 1
        F = (1 .- P) .* F_0 .+ P .* F_1
        EV = inv(I - β * F) * ((1 .- P).*(U_0 + e_0) + P .* (U_1 + e_1))
        V_tilde = (U_1 .+ β .* F_1 * EV) - (U_0 .+ β .* F_0 * EV)
        P_next = (1 .+ exp.(-V_tilde)).^(-1)
        err = maximum(abs.(P .- P_next))
        P = P_next
        # println("Iteration #", i)
        # println("Current error:", err)
    end

    simdata_probability = zeros(size(simdata)[1])
    simdata_choice = simdata.choice

    for i in 1:size(simdata)[1]
        simdata_probability[i] = P[simdata.state_id[i]+1]
    end

    sum(simdata_choice .* log.(simdata_probability) .+ (1 .- simdata_choice) .* log.(1 .- simdata_probability))
end

λ_grid = -8:0.1:0
log_likelihoods = zeros(length(λ_grid))

for i in 1:length(λ_grid)
    log_likelihoods[i] = log_likelihood(λ_grid[i], simdata, estimated_frequencies, state_space_df, F_0, F_1)
end

plot(λ_grid, log_likelihoods, legend = false)
ylabel!("Log-Likelihood")
xlabel!("λ")

savefig("question_4.png")

optimize(λ_temp -> -log_likelihood(λ_temp, simdata, estimated_frequencies, state_space_df, F_0, F_1), -10, 0, Brent())
