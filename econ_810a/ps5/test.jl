# Alex von Hafften
# ECON 810: Advanced Macro
# PS 5 - Code for testing the model
# Professor Carter Braxton

using Plots

cd("/Users/vonhafften/Documents/uw_madison/problem_sets/econ_810a/ps5/")

include("model.jl")

P = Primitives()
R = Initialize()

# test t = 12

solve_12!(R)
plot(R.vf_12, legend = false)

# test t = 11

solve_11!(R)
plot(R.vf_11, legend = false)
plot(R.pf_b_11, legend = false)

# test t = 10

solve_10!(R)
plot(R.vf_10)
plot(R.pf_b_10)

# test t = 9

solve_9!(R)
plot(R.vf_9[:, :, 1])
plot(R.vf_9[:, :, P.N_h]) # should be the same
plot(R.vf_9[:, 1, :])
plot(R.vf_9[:, P.N_h, :]) # should be the same
plot(R.pf_b_9[:, :, 1])
plot(R.pf_b_9[:, :, P.N_h]) # should be the same
plot(R.pf_τ_9[:, :, 1]) # should be zero because vf_4 is inititalized at zero
plot(R.pf_τ_9[:, :, P.N_h]) # should be the same

# test t = 8

solve_8!(R)
plot(R.vf_8[:, :, 1])
plot(R.vf_8[:, :, P.N_h]) # should be the same
plot(R.vf_8[:, 1, :])
plot(R.vf_8[:, P.N_h, :]) # should be the same
plot(R.pf_b_8[:, :, 1])
plot(R.pf_b_8[:, :, P.N_h]) # should be the same
plot(R.pf_i_8[:, :, 1]) # should be zero because vf_4 is inititalized at zero
plot(R.pf_i_8[:, :, P.N_h]) # should be the same
plot(R.pf_i_8[:, 1, :]) # should be zero because vf_4 is inititalized at zero
plot(R.pf_i_8[:, P.N_h, :]) # should be the same

# test t = 7

solve_7!(R)
plot(R.vf_7[:, :, 1])
plot(R.vf_7[:, :, P.N_h]) # should be the same
plot(R.vf_7[:, 1, :])
plot(R.vf_7[:, P.N_h, :]) # should be the same
plot(R.pf_b_7[:, :, 1])
plot(R.pf_b_7[:, :, P.N_h]) # should be the same
plot(R.pf_b_7[:, 1, :])
plot(R.pf_b_7[:, P.N_h, :]) # should be the same
plot(R.pf_i_7[:, :, 1]) # should be zero because vf_4 is inititalized at zero
plot(R.pf_i_7[:, :, P.N_h]) # should be the same

# test t = 6

solve_6!(R)
plot(R.vf_6[:, :, 1])
plot(R.vf_6[:, :, P.N_h]) # should be the same
plot(R.vf_6[:, 1, :])
plot(R.vf_6[:, P.N_h, :]) # should be the same
plot(R.pf_b_6[:, :, 1])
plot(R.pf_b_6[:, :, P.N_h]) # should be the same
plot(R.pf_b_6[:, 1, :])
plot(R.pf_b_6[:, P.N_h, :]) # should be the same
plot(R.pf_i_6[:, :, 1]) # should be zero because vf_4 is inititalized at zero
plot(R.pf_i_6[:, :, P.N_h]) # should be the same

# test t = 5

solve_5!(R)
plot(R.vf_5[:, :, 1])
plot(R.vf_5[:, :, P.N_h]) # should be the same
plot(R.vf_5[:, 1, :])
plot(R.vf_5[:, P.N_h, :]) # should be the same
plot(R.pf_b_5[:, :, 1])
plot(R.pf_b_5[:, :, P.N_h]) # should be the same
plot(R.pf_b_5[:, 1, :])
plot(R.pf_b_5[:, P.N_h, :]) # should be the same
plot(R.pf_i_5[:, :, 1]) # should be zero because vf_4 is inititalized at zero
plot(R.pf_i_5[:, :, P.N_h]) # should be the same

# test t = 4

solve_4!(R)
plot(R.vf_4[:, :])
plot(R.pf_b_4[:, :])

# test whole iterated problem

R = Initialize()
Solve!(R)

# test t = 12

plot(R.vf_12)

# test t = 11

plot(R.vf_11)
plot(R.pf_b_11)

# test t = 10

plot(R.vf_10)
plot(R.pf_b_10)

# test t = 9

plot(R.vf_9[:, :, 1])
plot(R.vf_9[:, :, P.N_h]) # should be the same
plot(R.vf_9[:, 1, :])
plot(R.vf_9[:, P.N_h, :]) # should be the same
plot(R.pf_b_9[:, :, 1])
plot(R.pf_b_9[:, :, P.N_h]) # should be the same
plot(R.pf_τ_9[:, :, 1]) # should be zero because vf_4 is inititalized at zero
plot(R.pf_τ_9[:, :, P.N_h]) # should be the same

# test t = 8

plot(R.vf_8[:, :, 1])
plot(R.vf_8[:, :, P.N_h]) # should be the same
plot(R.vf_8[:, 1, :])
plot(R.vf_8[:, P.N_h, :]) # should be the same
plot(R.pf_b_8[:, :, 1])
plot(R.pf_b_8[:, :, P.N_h]) # should be the same
plot(R.pf_i_8[:, :, 1]) # should be zero because vf_4 is inititalized at zero
plot(R.pf_i_8[:, :, P.N_h]) # should be the same
plot(R.pf_i_8[:, 1, :]) # should be zero because vf_4 is inititalized at zero
plot(R.pf_i_8[:, P.N_h, :]) # should be the same


# test t = 7

plot(R.vf_7[:, :, 1])
plot(R.vf_7[:, :, P.N_h]) # should be the same
plot(R.vf_7[:, 1, :])
plot(R.vf_7[:, P.N_h, :]) # should be the same
plot(R.pf_b_7[:, :, 1])
plot(R.pf_b_7[:, :, P.N_h]) # should be the same
plot(R.pf_b_7[:, 1, :])
plot(R.pf_b_7[:, P.N_h, :]) # should be the same
plot(R.pf_i_7[:, :, 1]) # should be zero because vf_4 is inititalized at zero
plot(R.pf_i_7[:, :, P.N_h]) # should be the same
plot(R.pf_i_7[:, 1, :]) # should be zero because vf_4 is inititalized at zero
plot(R.pf_i_7[:, P.N_h, :]) # should be the same

# test t = 6

plot(R.vf_6[:, :, 1])
plot(R.vf_6[:, :, P.N_h]) # should be the same
plot(R.vf_6[:, 1, :])
plot(R.vf_6[:, P.N_h, :]) # should be the same
plot(R.pf_b_6[:, :, 1])
plot(R.pf_b_6[:, :, P.N_h]) # should be the same
plot(R.pf_b_6[:, 1, :])
plot(R.pf_b_6[:, P.N_h, :]) # should be the same
plot(R.pf_i_6[:, :, 1]) # should be zero because vf_4 is inititalized at zero
plot(R.pf_i_6[:, :, P.N_h]) # should be the same

# test t = 5

plot(R.vf_5[:, :, 1])
plot(R.vf_5[:, :, P.N_h]) # should be the same
plot(R.vf_5[:, 1, :])
plot(R.vf_5[:, P.N_h, :]) # should be the same
plot(R.pf_b_5[:, :, 1])
plot(R.pf_b_5[:, :, P.N_h]) # should be the same
plot(R.pf_b_5[:, 1, :])
plot(R.pf_b_5[:, P.N_h, :]) # should be the same
plot(R.pf_i_5[:, :, 1]) # should be zero because vf_4 is inititalized at zero
plot(R.pf_i_5[:, :, P.N_h]) # should be the same

# test t = 4

plot(R.vf_4[:, :])
plot(R.pf_b_4[:, :])