# Code to solve Hennessy and Whited 2007

# Alex von Hafften
# Sept. 27, 2022

using Parameters, QuantEcon, Interpolations, Optim

include("structures.jl")
include("helper_functions.jl")

# applies bellman operator
function Apply_Bellman!(P::Primitives, G::Grids, R::Results)

    # interpolate interest rates and value
    q_interp = LinearInterpolation((G.grid_k, G.grid_b, G.grid_lz), R.q)
    vf_interp = LinearInterpolation((R.grid_w, G.grid_lz), R.vf)

    # initialize next guess for value function
    vf_next = zeros(R.N_w, G.N_lz)

    # iterate over state variables
    for (i_w, w) = enumerate(R.grid_w), (i_lz, lz) = enumerate(G.grid_lz)

        candidate_max = -Inf

        for (i_k_p, k_p) = enumerate(G.grid_k), (i_b_p, b_p) = enumerate(G.grid_b)
            
            # bond price
            q = q_interp(k_p, b_p, lz)

            candidate = w + q*b_p - k_p           # cash dividend if positive and equity issuance if negative
            candidate -= T_d(w + q*b_p - k_p, P)  # cash dividend tax
            candidate -= Λ(-(w + q*b_p - k_p), P) # equity issuance cost
 
            # add continuation value
            for (i_lz_p, lz_p) = enumerate(G.grid_lz)
                y_p = exp(lz_p) * k_p ^ P.α - P.δ * k_p - (1-q) * b_p     # taxable income
                w_p = y_p - T_C(y_p, P) + k_p - q*b_p                     # realized net worth
                candidate += 1/(1+P.r*(1 - P.τ_i)) * G.Π_lz[i_lz, i_lz_p] * max(vf_interp(w_p, lz_p), 0.0)
            end

            if (candidate > candidate_max)
                candidate_max = candidate
                vf_next[i_w, i_lz] = candidate
                R.pf_b[i_w, i_lz] = b_p
                R.pf_k[i_w, i_lz] = k_p
            end
        end
    end
    return vf_next
end


# solves for policy and value function with bond prices given
function Solve_vf!(P::Primitives, G::Grids, R::Results; show_progress::Bool = true)
    err = 100
    tolerence = 0.0001
    i = 1
    max_iterations = 10
    
    while true
        if show_progress
            println(err)
        end
        vf_next = Apply_Bellman!(P, G, R)
        err = maximum(abs.(vf_next - R.vf))
        R.vf = vf_next
        if err < tolerence
            break
        end
        if i > max_iterations
            println("Max iterations hit in Solve_vf!")
            break
        end
        i+=1
    end
end

function compute_default_states!(P::Primitives, G::Grids, R::Results)
    # determine default net worth threshold
    for i_lz = 1:G.N_lz
        # vf must be monotone
        if sum((R.vf[1:end-1, i_lz] - R.vf[2:end, i_lz]) .< 0.0) != R.N_w - 1
            println("R.vf is not monotone")
        end
        # invert value function
        vf_inv_interp = LinearInterpolation(R.vf[:, i_lz], R.grid_w)
        # get value at zero
        R.w_bar[i_lz] = vf_inv_interp(0.0)
    end

    # determine default threshold
    for (i_k_p, k_p) = enumerate(G.grid_k), (i_b_p, b_p) = enumerate(G.grid_b), (i_lz, lz) = enumerate(G.grid_lz)
        q = R.q[i_k_p, i_b_p, i_lz] # get bond price

        # compute net worth
        y_p = exp.(G.grid_lz) .* k_p .^ P.α .- P.δ .* k_p .- (1 .- q) .* b_p
        w_p = y_p .- broadcast(y_p -> T_C(y_p, P), y_p) .+ k_p .- q .* b_p

        # net worth over default threshold
        w_over_w_bar = w_p - R.w_bar

        # always under threshold set to max lz
        if sum(w_over_w_bar .< 0.0) == G.N_lz
            R.lz_d[i_k_p, i_b_p, i_lz] = G.max_lz
        elseif sum(w_over_w_bar .< 0.0) == 0
            R.lz_d[i_k_p, i_b_p, i_lz] = G.min_lz
        else
            # w_over_w_bar must be monotone
            if sum((w_over_w_bar[1:end-1] - w_over_w_bar[2:end]) .< 0.0) != G.N_lz - 1
                println("w_over_w_bar is not monotone")
            end
            # invert value function
            w_over_w_bar_inv_interp = LinearInterpolation(w_over_w_bar, G.grid_lz)
            # get value at zero
            R.lz_d[i_k_p, i_b_p, i_lz] = w_over_w_bar_inv_interp(0.0)
        end
    end
end

function compute_q!(P::Primitives, G::Grids, R::Results)
   
    q_next = zeros(G.N_k, G.N_b, G.N_lz)

    for (i_k_p, k_p) = enumerate(G.grid_k), (i_b_p, b_p) = enumerate(G.grid_b), (i_lz, lz) = enumerate(G.grid_lz)
        recovery_return = 0.0
        nondefault_mass = 0.0
        for (i_lz_p, lz_p) = enumerate(G.grid_lz)
            if lz_p >= R.lz_d[i_k_p, i_b_p, i_lz]
                nondefault_mass += G.Π_lz[i_lz, i_lz_p]
            else
                recovery_value = (1-P.ξ)*(1-P.δ)*k_p + exp(lz_p) * k_p^P.α - T_C(exp(lz_p) * k_p^P.α - P.δ * k_p, P) #- R.w_bar[i_lz_p]
                if recovery_value < 0.0
                    error("recovery value negative")
                end
                recovery_return += recovery_value / b_p * G.Π_lz[i_lz, i_lz_p]
            end
        end
        # compute interest rate
        r_tilde_next = (1/(1-P.τ_i)) * ((1 + P.r * (1 - P.τ_i) - recovery_return)/nondefault_mass - 1)

        # convert to bond price
        q_next[i_k_p, i_b_p, i_lz] = max(min(1/(1+r_tilde_next), 1.0), 0.0)
    end
    return q_next
end

function Solve_model(;show_progress = true)
    P = Initialize_Primitives()
    G = Initialize_Grids(P)
    R = Initialize_Results(P, G)
    
    tolerence = 0.001
    err = 100.0
    i = 1
    max_iterations = 200

    λ = 0.9

    while true
        Solve_vf!(P, G, R; show_progress = false)
        compute_default_states!(P, G, R)
        q_guess = compute_q!(P, G, R)
        q_next = R.q .* λ .+ q_guess .* (1-λ);
        err = maximum(abs.(R.q - q_next))
        R.q = q_next
        if show_progress
            println(err)
        end

        if err < tolerence
            break
        end
        if i > max_iterations
            println("Max iterations hit in Solve_model")
            break
        end

        # update net worth grid
        R.grid_w = compute_w_grid(R.q, R.N_w, P, G)
        R.min_w  = minimum(R.grid_w)
        R.max_w  = maximum(R.grid_w)

        i+=1
    end
    return R
end