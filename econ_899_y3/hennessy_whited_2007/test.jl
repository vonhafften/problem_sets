



########################################################
# test helper_functions
########################################################

using Plots

include("structures.jl");
include("helper_functions.jl");

P = Initialize_Primitives()
G = Initialize_Grids(P)

# test compute_T_C at r_tilde = r
grid_x = -100.0:1.0:100.0
grid_T_C = zeros(length(grid_x))
grid_T_d = zeros(length(grid_x))
for (i_x, x) = enumerate(grid_x)
    grid_T_C[i_x] = T_C(x, P)
    grid_T_d[i_x] = T_d(x, P)
end

plot(grid_x, grid_T_C) 
plot(grid_x, grid_T_d) 

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

# find defaulting states
R.pf_d[R.vf .<= 0.0] .= 1

# iterate over productiviies
for i_lz = 1:G.N_lz_c
    # vf must be 
    if sum((R.vf[1:end-1, i_lz] - R.vf[2:end, i_lz]) .< 0.0) != G.N_w - 1
        println("R.vf is not monotone")
    end

    vf_inv_interp = LinearInterpolation(R.vf[:, i_lz], G.grid_w)
    R.w_bar[i_lz] = vf_inv_interp(0.0)
end

for (i_k_p, k_p) = enumerate(G.grid_k), (i_b_p, b_p) = enumerate(G.grid_b), (i_lz, lz) = enumerate(G.grid_lz)
    R.r_tilde[i_k_p, i_b_p, i_lz] = (1/(1-P.τ_i)) * ((1 + P.r*(1-P.τ_i) - )/())
end
