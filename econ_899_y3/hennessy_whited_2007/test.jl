cd(@__DIR__)

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

plot(R.grid_w, R.vf)
plot(R.grid_w, R.pf_b)
plot(R.grid_w, R.pf_k)
plot!(R.grid_w, R.grid_w)

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
compute_default_states!(P, G, R)
q_next = compute_q!(P, G, R)
err = sum(abs.(R.q - q_next))
R.q = 0.9 .* R.q .+ 0.1 .* q_next

R.grid_w = compute_w_grid(R.q, R.N_w, P, G)
R.min_w  = minimum(R.grid_w)
R.max_w  = maximum(R.grid_w)

plot(G.grid_k, R.q[:,G.N_b,:])



Solve_vf!(P, G, R)
compute_default_states!(P, G, R)
q_next = compute_q!(P, G, R)
err = sum(abs.(R.q - q_next))
R.q = q_next;

Solve_vf!(P, G, R2)
compute_default_states!(P, G, R2)
q_next = compute_q!(P, G, R2)
err = sum(abs.(R2.q - q_next))
R2.q = q_next;

R1.q .-R2.q
plot(G.grid_lz, R.w_bar)
plot!(G.grid_lz, R2.w_bar)


plot(R1.grid_w, R.vf)
plot!(R2.grid_w, R2.vf)

plot(G.grid_k, R.lz_d[:,G.N_b,:])
plot!(G.grid_k, R2.lz_d[:,G.N_b,:])

plot(G.grid_k, R2.q[:,G.N_b,:])



plot(R.grid_w[1:50], R1.vf[1:50, :])
plot!(R.grid_w[1:50], R2.vf[1:50, :])
plot(R.grid_w, R.pf_b)
plot(R.grid_w, R.pf_k)   # on second pass this is weird. Low productivity firms choose more capital.
plot(G.grid_lz, R.pf_k') # on second pass this is weird. Low productivity firms choose more capital.
plot(G.grid_lz, R.w_bar) # on second pass this is weird. Low productivity firms choose more capital.

plot(G.grid_k, R.lz_d[:,G.N_b,:])
plot(G.grid_k, R.lz_d[:,G.N_b-1,:])
plot(G.grid_k, R.q[:,G.N_b,:])
plot(G.grid_k, R.q[:,G.N_b-1,:])


Solve_vf!(P, G, R)
compute_default_states!(P, G, R)
q_next = compute_q!(P, G, R)
err = sum(abs.(R.q - q_next))
R.q = q_next;

plot(R.grid_w, R.vf)
plot(R.grid_w[1:50], R.vf[1:50, :])
plot(R.grid_w, R.pf_b)
plot(G.grid_lz, R.pf_b') # on second pass this is weird. Low productivity firms choose more capital.
plot(R.grid_w, R.pf_k) # on second pass this is weird. Low productivity firms choose more capital.
plot(G.grid_lz, R.pf_k') # on second pass this is weird. Low productivity firms choose more capital.
plot(G.grid_lz, R.w_bar) # on second pass this is weird. Low productivity firms choose more capital.

plot(G.grid_k, R.lz_d[:,G.N_b,:])
plot(G.grid_k, R.lz_d[:,G.N_b-1,:])
plot(G.grid_k, R.lz_d[:,G.N_b-2,:])
plot(G.grid_k, R.q[:,G.N_b,:])
plot(G.grid_k, R.q[:,G.N_b-1,:])
plot(G.grid_k, R.q[:,G.N_b-2,:])
plot(G.grid_k, R.q[:,G.N_b-3,:])
plot(G.grid_k, R.q[:,G.N_b-4,:])
plot(G.grid_k, R.q[:,G.N_b-5,:])
plot(G.grid_k, R.q[:,G.N_b-6,:])


########################################################
# test model.jl
########################################################

using Plots

include("structures.jl");
include("helper_functions.jl");
include("model.jl");

R = Solve_model()

plot(R.grid_w, R.vf)
plot(R.grid_w, R.pf_b)
plot(G.grid_lz, R.pf_b') # on second pass this is weird. Low productivity firms choose more capital.
plot(R.grid_w, R.pf_k) # on second pass this is weird. Low productivity firms choose more capital.
plot(G.grid_lz, R.pf_k') # on second pass this is weird. Low productivity firms choose more capital.
plot(G.grid_lz, R.w_bar) # on second pass this is weird. Low productivity firms choose more capital.

plot(G.grid_k, R.lz_d[:,G.N_b,:])
plot(G.grid_k, R.lz_d[:,G.N_b-1,:])
plot(G.grid_k, R.lz_d[:,G.N_b-2,:])
plot(G.grid_k, R.lz_d[:,G.N_b-3,:])
plot(G.grid_k, R.lz_d[:,G.N_b-4,:])
plot(G.grid_k, R.lz_d[:,G.N_b-5,:])
plot(G.grid_k, R.lz_d[:,G.N_b-6,:])

plot(G.grid_k, R.q[:,G.N_b,:])
plot(G.grid_k, R.q[:,G.N_b-1,:])
plot(G.grid_k, R.q[:,G.N_b-2,:])
plot(G.grid_k, R.q[:,G.N_b-3,:])
plot(G.grid_k, R.q[:,G.N_b-4,:])
plot(G.grid_k, R.q[:,G.N_b-5,:])
plot(G.grid_k, R.q[:,G.N_b-6,:])