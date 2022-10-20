
include("model.jl")

# Initialize objects
P = Primitives()
G = Grids()
R = Initialize_Results(P, G)

# solve vf for initial bond price guess
solve_vf!(P, G, R)

# figures
plot(G.grid_a, R.vf_h0', label = "y = " .* string.(G.grid_y') .* "; h = 0", title = "VF", color = ["red" "blue"]);
plot!(G.grid_a_p, R.vf_h1', label = "y = " .* string.(G.grid_y') .* "; h = 1", color = ["red" "blue"], linestyle = :dash)

plot(G.grid_a, R.pf_d_h0', label = "y = " .* string.(G.grid_y') .* "; h = 0", title = "PF d", color = ["red" "blue"])

plot(G.grid_a, R.pf_a_h0', label = "y = " .* string.(G.grid_y') .* "; h = 0", title = "PF a", color = ["red" "blue"]);
plot!(G.grid_a_p, R.pf_a_h1', label = "y = " .* string.(G.grid_y') .* "; h = 1", color = ["red" "blue"], linestyle = :dash);
plot!(G.grid_a, G.grid_a, label = "45ᵒ line", color = "black", legend = :bottomright)

plot(G.grid_a_n, R.q', label = "y = " .* string.(G.grid_y'), title = "q", color = ["red" "blue"])

# test computing q
q_next = compute_q(P, G, R)
plot(G.grid_a_n, q_next', label = "y = " .* string.(G.grid_y'), title = "q", color = ["red" "blue"])

# test solving for q
solve_q!(P, G, R)

plot(G.grid_a, R.vf_h0', label = "y = " .* string.(G.grid_y') .* "; h = 0", title = "VF", color = ["red" "blue"]);
plot!(G.grid_a_p, R.vf_h1', label = "y = " .* string.(G.grid_y') .* "; h = 1", color = ["red" "blue"], linestyle = :dash, legend = :bottomright)
savefig("vf.png")

plot(G.grid_a, R.pf_d_h0', label = "y = " .* string.(G.grid_y') .* "; h = 0", title = "PF d", color = ["red" "blue"])
savefig("pf_d.png")

plot(G.grid_a, R.pf_a_h0', label = "y = " .* string.(G.grid_y') .* "; h = 0", title = "PF a", color = ["red" "blue"]);
plot!(G.grid_a_p, R.pf_a_h1', label = "y = " .* string.(G.grid_y') .* "; h = 1", color = ["red" "blue"], linestyle = :dash);
plot!(G.grid_a, G.grid_a, label = "45ᵒ line", color = "black", legend = :bottomright)
savefig("pf_a.png")

plot(G.grid_a_n, R.q', label = "y = " .* string.(G.grid_y'), title = "q", color = ["red" "blue"], legend = :topleft)
savefig("q.png")
