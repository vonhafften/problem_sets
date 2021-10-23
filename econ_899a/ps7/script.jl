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

function J_TH(b, which_moments, W)

    ρ = b[1] # unpack b
    σ = b[2]

    x = zeros(T, H) # initialize x matrix

    # fill in x based on parameters
    for i_H = 1:H
        x[1, i_H] = e[1, i_H]
        for t = 2:T
            x[t, i_H] = ρ * x[t - 1, i_H] + σ * e[t, i_H]
        end
    end

    x_bar = repeat(sum(x; dims = 1)./T, T) # compute mean by simulation
    x_lag = vcat(zeros(H)', x[1:T-1, :])   # lag x

    m = zeros(T, H, 3)
    m[:,:,1] = x
    m[:,:,2] = (x .- x_bar).^2
    m[:,:,3] = (x .- x_bar) .* (x_lag .- x_bar)
    
    # compute M_TH
    M_TH = reshape(sum(m; dims = 1:2)/(T*H), 3, 1)

    # evaluate J_TH
    (M_T[which_moments] - M_TH[which_moments])' * W * (M_T[which_moments] - M_TH[which_moments])
end

# plotting function
function J_TH_plotting(ρ, σ)
    return J_TH([ρ, σ], 1:2, I)
end

# grids
ρ_grid = 0.35:0.01:0.65;
σ_grid = 0.8:0.01:1.2;

# Surface of J_TH plot
surface(ρ_grid, σ_grid, J_TH_plotting);
xlabel!("ρ");
ylabel!("σ")
savefig("figures/q4a.png");

# b_TH_hat_1
opt = optimize(b -> J_TH(b, 1:2, I), [ρ_0, σ_0]);
ρ_hat_1 = Optim.minimizer(opt)[1]
σ_hat_1 = Optim.minimizer(opt)[2]

##########################################################################################
# Question 4b
##########################################################################################

function compute_W_hat_star(which_moments)
    n_moments = length(which_moments)

    # simulate x_1
    x_1 = zeros(T, H) 
    for i_H = 1:H
        x_1[1, i_H] = e[1, i_H]
        for t = 2:T
            x_1[t, i_H] = ρ_hat_1 * x_1[t - 1, i_H] + σ_hat_1 * e[t, i_H]
        end
    end

    x_1_bar = repeat(sum(x_1; dims = 1)./T, T) # compute mean by simulation
    x_1_lag = vcat(zeros(H)', x_1[1:T-1, :])   # lag x

    m_1 = zeros(T, H, 3)
    m_1[:,:,1] = x_1
    m_1[:,:,2] = (x_1 .- x_1_bar).^2
    m_1[:,:,3] = (x_1 .- x_1_bar) .* (x_1_lag .- x_1_bar)

    # compute M_TH
    M_TH_1 = reshape(sum(m_1; dims = 1:2)/(T*H), 3)

    function compute_Γ(j)
        Γ = zeros(3, 3)

        for i_1 = 1:3
            for i_2 in 1:3
                Γ[i_1, i_2] = sum((m_1[j+1:T,:,i_1] .- M_TH_1[i_1]) * (m_1[1:T-j, :, i_2] .- M_TH_1[i_2])')/ (T * H)
            end
        end
        Γ[which_moments, which_moments]
    end

    S_TH = compute_Γ(0)

    lags = 4

    for j = 1:lags
        Γ_j = compute_Γ(j)
        S_TH += (1 - (j/(lags + 1))) * (Γ_j + Γ_j')
    end

    S_TH = (1 + 1/H) * S_TH

    inv(S_TH)
end

W_hat_star = compute_W_hat_star(1:2)

opt = optimize(b -> J_TH(b, 1:2, W_hat_star), [ρ_0, σ_0]);
ρ_hat_2 = Optim.minimizer(opt)[1]
σ_hat_2 = Optim.minimizer(opt)[2]

##########################################################################################
# Question 4c
##########################################################################################



##########################################################################################
# Question 4d
##########################################################################################


##########################################################################################
# Question 5a
##########################################################################################

function J_TH_plotting(ρ, σ)
    return J_TH([ρ, σ], 2:3, I)
end

# Surface of J_TH plot
surface(ρ_grid, σ_grid, J_TH_plotting);
xlabel!("ρ");
ylabel!("σ")
savefig("figures/q5a.png");

# b_hat_1
opt = optimize(b -> J_TH(b, 2:3, I), [ρ_0, σ_0]);
ρ_hat_1 = Optim.minimizer(opt)[1]
σ_hat_1 = Optim.minimizer(opt)[2]

##########################################################################################
# Question 5b
##########################################################################################

W_hat_star = compute_W_hat_star(2:3)

opt = optimize(b -> J_TH(b, 1:2, W_hat_star), [ρ_0, σ_0]);
ρ_hat_2 = Optim.minimizer(opt)[1]
σ_hat_2 = Optim.minimizer(opt)[2]

##########################################################################################
# Question 5c
##########################################################################################



##########################################################################################
# Question 5d
##########################################################################################



##########################################################################################
# Question 5a
##########################################################################################

function J_TH_plotting(ρ, σ)
    return J_TH([ρ, σ], 1:3, I)
end

# Surface of J_TH plot
surface(ρ_grid, σ_grid, J_TH_plotting_23);
xlabel!("ρ");
ylabel!("σ")
savefig("figures/q5a.png");

# b_hat_1
opt = optimize(b -> J_TH(b, 1:3, I), [ρ_0, σ_0]);
ρ_hat_1 = Optim.minimizer(opt)[1]
σ_hat_1 = Optim.minimizer(opt)[2]

##########################################################################################
# Question 5b
##########################################################################################

W_hat_star = compute_W_hat_star(1:3)

opt = optimize(b -> J_TH(b, 1:3, W_hat_star), [ρ_0, σ_0]);
ρ_hat_2 = Optim.minimizer(opt)[1]
σ_hat_2 = Optim.minimizer(opt)[2]

##########################################################################################
# Question 5c
##########################################################################################



##########################################################################################
# Question 5d
##########################################################################################

