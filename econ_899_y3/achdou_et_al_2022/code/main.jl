# Code to solve Huggett model in cts time 
# a la Achdou et al (2017)

# Alex von Hafften
# Sept. 22

using Plots

cd(@__DIR__)

include("model.jl")

R_low_r = Initialize(0.005, -0.3)
Solve!(R_low_r)

R_high_r = Initialize(0.01, -0.3)
Solve!(R_high_r)

plot(R_low_r.a_grid, R_low_r.vf, label = ["low r and low y" "low r and high y"], color = :blue, line = [:solid :dash])
plot!(R_high_r.a_grid, R_high_r.vf, label = ["high r and low y" "high r and high y"], color = :red, line = [:solid :dash])

plot(R_low_r.a_grid, R_low_r.π, label = ["low r and low y" "low r and high y"], color = :blue, line = [:solid :dash])
plot!(R_high_r.a_grid, R_high_r.π, label = ["high r and low y" "high r and high y"], color = :red, line = [:solid :dash])

plot(R_low_r.a_dot, label = ["low r and low y" "low r and high y"], color = :blue, line = [:solid :dash])
plot!(R_high_r.a_dot, label = ["high r and low y" "high r and high y"], color = :red, line = [:solid :dash])

# solve for equilibrium interest rate with market clearing
R_loose_bc = Solve!(0.0001, 0.01, -0.3);
R_tight_bc = Solve!(0.0001, 0.01, -0.2);

plot(R_loose_bc.a_grid, R_loose_bc.vf, label = ["loose bc; low y" "loose bc; high y"], color = :blue, line = [:solid :dash])
plot!(R_tight_bc.a_grid, R_tight_bc.vf, label = ["tight bc; low y" "tight bc; high y"], color = :red, line = [:solid :dash])

plot(R_loose_bc.a_grid, R_loose_bc.π, label =  ["loose bc; low y" "loose bc; high y"], color = :blue, line = [:solid :dash])
plot!(R_tight_bc.a_grid, R_tight_bc.π, label = ["tight bc; low y" "tight bc; high y"], color = :red, line = [:solid :dash])

plot(R_loose_bc.a_grid, R_loose_bc.a_dot, label = ["loose bc; low y" "loose bc; high y"], color = :blue, line = [:solid :dash])
plot!(R_tight_bc.a_grid, R_tight_bc.a_dot, label = ["tight bc; low y" "tight bc; high y"], color = :red, line = [:solid :dash])


