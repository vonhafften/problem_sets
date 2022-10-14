
        else # solve using function minimizer
            # obj takes any real numbers as args and forces args to be in grid
            # return negative payoff (i.e. optimize finds minimum)
            function obj(args)
                # println(args)
                k_p = (G.max_k - G.min_k) * exp(args[1])/(1+exp(args[1])) + G.min_k
                b_p = (G.max_b - G.min_b) * exp(args[2])/(1+exp(args[2])) + G.min_b
                if (isnan(k_p))
                    k_p = G.max_k
                end
                if (isinf(k_p))
                    k_p = G.max_k
                end
                if (isnan(b_p))
                    b_p = G.max_b
                end
                if (isinf(b_p))
                    b_p = G.max_b
                end
                -payoff(k_p, b_p)
            end

            optim_obj = optimize(obj, [(G.max_k + G.min_k)/2, (G.max_b + G.min_b)/2], NelderMead())

            vf_next[i_w, i_lz] = -Optim.minimum(optim_obj)
            R.pf_k[i_w, i_lz]  = (G.max_k - G.min_k)* exp(optim_obj.minimizer[1])/(1+exp(optim_obj.minimizer[1])) +  G.min_k
            R.pf_b[i_w, i_lz]  = (G.max_b - G.min_b)* exp(optim_obj.minimizer[2])/(1+exp(optim_obj.minimizer[2])) +  G.min_b
        end



# find defaulting states
R.pf_d[R.vf .<= 0.0] .= 1

# determine default net worth threshold
# iterate over productivities
for i_lz = 1:G.N_lz_c
    # vf must be monotone
    if sum((R.vf[1:end-1, i_lz] - R.vf[2:end, i_lz]) .< 0.0) != R.N_w - 1
        println("R.vf is not monotone")
    end
    # invert value function
    vf_inv_interp = LinearInterpolation(R.vf[:, i_lz], R.grid_w)
    # get value at zero
    R.w_bar[i_lz] = vf_inv_interp(0.0)
end


r_tilde_next = zeros(G.N_k, G.N_b, G.N_lz_c)

for (i_k_p, k_p) = enumerate(G.grid_k), (i_b_p, b_p) = enumerate(G.grid_b), (i_lz, lz) = enumerate(G.grid_lz_c)
    payoff_f = 0.0
    payoff_s = 0.0

    r_tilde = R.r_tilde[i_k_p, i_b_p, i_lz]

    for (i_lz_p, lz_p) = enumerate(G.grid_lz_c)
        y_p = exp(lz_p) * k_p ^ P.α - P.δ * k_p - r_tilde * b_p # taxable income
        w_p = y_p - T_C(y_p, P) + k_p - b_p                     # realized net worth
        if w_p > R.w_bar[i_lz_p]
            payoff_s += G.Π_lz[i_lz, i_lz_p]
        else
            y_p = exp(lz_p) * k_p^P.α
            recovery_value = (1-P.ξ)*(1-P.δ)*k_p + y_p - T_C(y_p, P) -  R.w_bar[i_lz_p]
            payoff_f += recovery_value / b_p * G.Π_lz_c[i_lz, i_lz_p]
        end
    end

    r_tilde_next[i_k_p, i_b_p, i_lz] = (1/(1-P.τ_i)) * ((1 + P.r*(1-P.τ_i) - payoff_f)/payoff_s - 1)
end

# set artifical bound on bond rates
r_tilde_next[isinf.(r_tilde_next)] .= 1.0
r_tilde_next[r_tilde_next .> 1.0] .= 1.0

R.r_tilde = r_tilde_next

plot(R.r_tilde[:,:,10])

temp_grid_w = zeros(G.N_k, G.N_b, G.N_lz_c, G.N_lz_c)
for (i_k_p, k_p) = enumerate(G.grid_k), (i_b_p, b_p) = enumerate(G.grid_b), (i_lz, lz) = enumerate(G.grid_lz_c), (i_lz_p, lz_p) = enumerate(G.grid_lz_c)
    y_p = exp(lz_p) * k_p ^ P.α - P.δ * k_p - R.r_tilde[i_k_p,i_b_p,i_lz] * b_p
    w_p = y_p - T_C(y_p, P) + k_p - b_p
    temp_grid_w[i_k_p, i_b_p, i_lz, i_lz_p] = w_p
end
R.min_w  = minimum(temp_grid_w)
R.max_w  = maximum(temp_grid_w)
R.grid_w = collect(range(R.min_w, R.max_w; length = R.N_w))

Solve_vf!(P, G, R)

plot(R.grid_w, R.vf)
plot(R.grid_w, R.pf_b)
plot(R.grid_w, R.pf_k)s