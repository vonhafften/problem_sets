##########################################################################################
# Problem set 7
# ECON 899A: Computational Economics
#
# Alex von Hafften
# October 27, 2021
##########################################################################################

using Random, Distributions, Plots, LinearAlgebra, Interpolations, Optim

# Parameters
ρ_0 = 0.5;
σ_0 = 1.0;
T   = 200;
H   = 10;

dist = Normal(0, 1.0);

##########################################################################################
# Question 2
##########################################################################################

# Draw shocks
Random.seed!(1234);
ε = rand(dist, T) * σ_0;

# Create data
x_0 = zeros(T);
x_0[1] = ε[1];
for t = 2:T
    x_0[t] = ρ_0 * x_0[t - 1] + ε[t]
end

plot(x_0, legend = false)
savefig("figures/q2.png");

# compute M_T
M_T = zeros(3);
M_T[1] = sum(x_0) / T
M_T[2] = sum((x_0 .- M_T[1]).^2) / T
M_T[3] = sum((x_0 .- M_T[1]).* (prepend!(x_0[1:T-1], 0) .- M_T[1])) / T

##########################################################################################
# Question 3
##########################################################################################

Random.seed!(2345);
e = reshape(rand(dist, T * H), T, H);

##########################################################################################
# Question 4a
##########################################################################################

# grids
ρ_grid = 0.35:0.01:0.65;
σ_grid = 0.8:0.01:1.2;

# initialize vectors
x = zeros(T, H);
M_TH = zeros(2);
J_TH = zeros(length(ρ_grid), length(σ_grid));

# iterate over each ρ and σ
for (i_ρ, ρ) = enumerate(ρ_grid)
    for (i_σ, σ) = enumerate(σ_grid)

        # loop over each trial
        for i_H = 1:H
            x[1, i_H] = e[1, i_H]
            for t = 2:T
                x[t, i_H] = ρ * x[t - 1, i_H] + σ * e[t, i_H]
            end
        end

        # compute M_TH
        M_TH[1] = sum(x) / (T*H)
        M_TH[2] = sum((x .- M_TH[1]).^2) / (T*H)

        # evaluate J_TH
        J_TH[i_ρ, i_σ] = (M_T[1:2] - M_TH)' * I * (M_T[1:2] - M_TH)
    end
end

# Cubic interpolate between points; scale onto grids
J_TH_interp = interpolate(J_TH, BSpline(Cubic(Line(OnGrid()))));
J_TH_interp = Interpolations.scale(J_TH_interp, ρ_grid,  σ_grid);

# a function to plot the surface so that the x and y coordinates are correct.
function J_TH_(b)
    ρ = b[1]
    σ = b[2]
    if ρ < minimum(ρ_grid) || ρ > maximum(ρ_grid) || σ < minimum(σ_grid) || σ > maximum(σ_grid)
        return 1/eps() # return a very large number if outside of grid
    end
    J_TH_interp(ρ, σ)
end

# specific plotting function
function J_TH_plotting(ρ, σ)
    return J_TH_([ρ, σ])
end

# Surface of J_TH plot
surface(ρ_grid, σ_grid, J_TH_plotting);
xlabel!("ρ");
ylabel!("σ")
savefig("figures/q4a.png");

# b_TH_hat_1
@elapsed opt = optimize(J_TH_, [ρ_0, σ_0])
ρ_hat_1 = Optim.minimizer(opt)[1]
σ_hat_1 = Optim.minimizer(opt)[2]

##########################################################################################
# Question 4b
##########################################################################################











##########################################################################################
# Question 5a
##########################################################################################

# grids
ρ_grid = 0.35:0.01:0.65;
σ_grid = 0.8:0.01:1.2;

# initialize vectors
x = zeros(T, H);
M_TH = zeros(3);
J_TH = zeros(length(ρ_grid), length(σ_grid));

# iterate over each ρ and σ
for (i_ρ, ρ) = enumerate(ρ_grid)
    for (i_σ, σ) = enumerate(σ_grid)

        # loop over each trial
        for i_H = 1:H
            x[1, i_H] = e[1, i_H]
            for t = 2:T
                x[t, i_H] = ρ * x[t - 1, i_H] + σ * e[t, i_H]
            end
        end

        # compute M_TH
        M_TH[1] = sum(x) / (T*H)
        M_TH[2] = sum((x .- M_TH[1]).^2) / (T*H)
        M_TH[3] = sum((x .- M_TH[1]).^2) / (T*H)  ##### Need to change

        # evaluate J_TH
        J_TH[i_ρ, i_σ] = (M_T[2:3] - M_TH[2:3])' * I * (M_T[2:3] - M_TH[2:3])
    end
end

# Cubic interpolate between points; scale onto grids
J_TH_interp = interpolate(J_TH, BSpline(Cubic(Line(OnGrid()))));
J_TH_interp = Interpolations.scale(J_TH_interp, ρ_grid,  σ_grid);

# a function to plot the surface so that the x and y coordinates are correct.
function J_TH_(b)
    ρ = b[1]
    σ = b[2]
    if ρ < minimum(ρ_grid) || ρ > maximum(ρ_grid) || σ < minimum(σ_grid) || σ > maximum(σ_grid)
        return 1/eps() # return a very large number if outside of grid
    end
    J_TH_interp(ρ, σ)
end

# specific plotting function
function J_TH_plotting(ρ, σ)
    return J_TH_([ρ, σ])
end

# Surface of J_TH plot
surface(ρ_grid, σ_grid, J_TH_plotting);
xlabel!("ρ");
ylabel!("σ")
savefig("figures/q4a.png");

# b_TH_hat_1
@elapsed opt = optimize(J_TH_, [ρ_0, σ_0])
ρ_hat_1 = Optim.minimizer(opt)[1]
σ_hat_1 = Optim.minimizer(opt)[2]