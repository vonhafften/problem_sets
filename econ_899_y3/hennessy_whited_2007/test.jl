



########################################################
# test helper_functions
########################################################

include("structures.jl");
include("helper_functions.jl");

P = Initialize_Primitives()
G = Initialize_Grids(P)

# test compute_π

π_grid = zeros(G.N_k)
for i = 1:G.N_k
    π_grid[i] = compute_π(G.grid_k[i], P)
end
plot(G.grid_k, π_grid)

# test compute_π_d_inv

π_d_inv_grid = zeros(G.N_k)
for i = 1:G.N_k
    π_d_inv_grid[i] = compute_π_d_inv(π_grid[i], P)
end
plot(G.grid_k, π_d_inv_grid)

# test compute_T_C at r_tilde = r

T_C_grid = zeros(G.N_k, G.N_b, G.N_lz_c)
for (i_k, k) = enumerate(G.grid_k), (i_b, b) = enumerate(G.grid_b), (i_lz, lz) = enumerate(G.grid_lz_c)
    T_C_grid[i_k, i_b, i_lz] = compute_T_C(k, b, exp(lz), P.r, P)
end

plot(G.grid_k, T_C_grid[:, [1, G.N_b],  Int64(floor(G.N_lz_c/2))]) # low and high amount of debt
plot(G.grid_k, T_C_grid[:, Int64(floor(G.N_b/2)), [1, G.N_lz_c]]) # low and high productivities
plot(G.grid_b, T_C_grid[[1, G.N_k], :, Int64(floor(G.N_lz_c/2))]') # low and high productivities

# test compute_T_d
x_grid = -10.0:1.0:100.0
T_d_grid = zeros(length(x_grid))
for (i_x, x) = enumerate(x_grid)
    T_d_grid[i_x] = compute_T_d(x, P)
end
plot(x_grid, T_d_grid)



########################################################
# test model.jl
########################################################

using Plots

include("structures.jl");
include("helper_functions.jl");
include("model.jl");

P = Initialize_Primitives()
G = Initialize_Grids(P)
R = Initialize_Results(P, G)
Solve_vf!(P, G, R)

plot(G.grid_w, R.vf)
plot(G.grid_w, R.pf_b)
plot(G.grid_w, R.pf_k)

# update defaulting states
R.pf_d[R.vf .<= 0.0] .= 1
for i_lz = 1:G.N_lz_c
    R.w_bar[i_lz] = maximum(G.grid_w[R.pf_d[:,i_lz] .== 1])
end

for (i_k_p, k_p) = enumerate(G.grid_k), (i_b_p, b_p) = enumerate(G.grid_b), (i_lz, lz) = enumerate(G.grid_lz)
    R.r_tilde[i_k_p, i_b_p, i_lz] = (1/(1-P.τ_i)) * ((1 + P.r*(1-P.τ_i) - )/())
end

compute_T_C(10.0, 10.0, 0.5, 0.01, P)





grid_2 = reshape(grid_w, N_k*N_b*N_lz_c*N_lz_c)



payoffs = zeros(G.N_k, G.N_b)
for (i_k, k) = enumerate(G.grid_k), (i_b, b) = enumerate(G.grid_b)
    payoffs[i_k, i_b] = payoff(k, b)
end