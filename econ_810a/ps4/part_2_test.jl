# Alex von Hafften
# ECON 810: Advanced Macro
# PS 4 - Part 2 - Code for testing the model code
# Professor Carter Braxton

using Plots, Statistics

cd("/Users/vonhafften/Documents/UW Madison/problem_sets/econ_810a/ps4/")

include("part_2_model.jl")

P = Primitives()
G = Grids()
R = Initialize(-0.029, 0.11, 2.0, 0.5)

Solve_terminal_period!(R)
plot(R.vf[P.T, [10, 50, 90], :]')
plot(R.vf[P.T, :, [10, 50, 90]])

Solve_nonterminal_period!(R)

# t = 29
which_t = 5

plot(R.vf[which_t, [10, 50, 90], :]')
plot(R.vf[which_t, :, [10, 50, 90]])


include("part_2_simulation.jl")

S = Initialize_simulation(1)
Simulate_model!(S, R)

plot(S.h')
plot(S.z')
plot(S.s')

S = Initialize_simulation(100000)
Simulate_model!(S, R)

# part a

plot(mean(S.k; dims= 1)')
plot(std(S.k; dims= 1)')

k_kurtosis = zeros(P.T)
for t in 1:P.T 
    k_kurtosis[t] = kurtosis(S.k[:, t])
end
plot(k_kurtosis)

# part b

which_t = 10

plot(R.pf_k[which_t, [10, 50, 90], :]')
plot(R.pf_k[which_t, :, [10, 50, 90]])

plot(R.pf_s[which_t, [10, 50, 90], :]')
plot(R.pf_s[which_t, :, [10, 50, 90]])

which_t = 20

plot(R.pf_k[which_t, [10, 50, 90], :]')
plot(R.pf_k[which_t, :, [10, 50, 90]])

plot(R.pf_s[which_t, [10, 50, 90], :]')
plot(R.pf_s[which_t, :, [10, 50, 90]])

# part c

R_2 = Initialize(-0.029, 0.22, 2.0, 0.5)
Solve!(R_2)

which_t = 10

plot(R_2.pf_k[which_t, [10, 50, 90], :]')
plot(R_2.pf_k[which_t, :, [10, 50, 90]])

plot(R_2.pf_s[which_t, [10, 50, 90], :]')
plot(R_2.pf_s[which_t, :, [10, 50, 90]])

which_t = 20

plot(R_2.pf_k[which_t, [10, 50, 90], :]')
plot(R_2.pf_k[which_t, :, [10, 50, 90]])

plot(R_2.pf_s[which_t, [10, 50, 90], :]')
plot(R_2.pf_s[which_t, :, [10, 50, 90]])

# part d

R_3 = Initialize(-0.029, 0.11, 2.0, 1.0)
