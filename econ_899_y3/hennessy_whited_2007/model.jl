# Code to solve Hennessy and Whited 2007

# Alex von Hafften
# Sept. 27, 2022

using Parameters, QuantEcon, Interpolations, Optim

include("structures.jl")
include("helper_functions.jl")

# applies bellman operator
function apply_Bellman!(P::Primitives, G::Grids, R::Results)

    # interpolate interest rates and value
    vf_interp = LinearInterpolation((G.grid_w, G.grid_lz), R.vf)

    # initialize next guess for value function
    vf_next = zeros(G.N_w, G.N_lz)

    # iterate over state variables
    for (i_w, w) = enumerate(G.grid_w), (i_lz, lz) = enumerate(G.grid_lz)

        candidate_max = -Inf

        for (i_k_p, k_p) = enumerate(G.grid_k), (i_b_p, b_p) = enumerate(G.grid_b)
            
            # bond price
            q = R.q[i_k_p, i_b_p, i_lz]

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
function solve_vf!(P::Primitives, G::Grids, R::Results; show_progress::Bool = true)
    err = 100
    tolerence = 0.0001
    i = 1
    max_iterations = 20
    
    while true
        if show_progress
            println(err)
        end
        vf_next = apply_Bellman!(P, G, R)
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
        if sum((R.vf[1:end-1, i_lz] - R.vf[2:end, i_lz]) .<= 0.0) != G.N_w - 1
            error("R.vf is not monotone")
        end

        # get value at zero by inverting  value function w linear Interpolations 
        # vf_inv_interp = LinearInterpolation(R.vf[:, i_lz], G.grid_w)
        # R.w_bar[i_lz] = vf_inv_interp(0.0)

        # use math; speeds up because you don't need the full interpolated object
        for i_w = 1:(G.N_w-1)
            if (R.vf[i_w, i_lz] < 0.0) & (R.vf[i_w+1, i_lz] > 0.0)
                m = (R.vf[i_w+1, i_lz] - R.vf[i_w, i_lz])/(G.grid_w[i_w+1] - G.grid_w[i_w])
                b = R.vf[i_w, i_lz] - m * G.grid_w[i_w]
                R.w_bar[i_lz] = -b/m
                break
            end
        end
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
            if sum((w_over_w_bar[1:end-1] - w_over_w_bar[2:end]) .<= 0.0) != G.N_lz - 1
                println(w_over_w_bar)
                error("w_over_w_bar is not monotone")
            end

            # invert value function to get value at zero using linear interpolation
            # w_over_w_bar_inv_interp = LinearInterpolation(w_over_w_bar, G.grid_lz)
            # R.lz_d[i_k_p, i_b_p, i_lz] = w_over_w_bar_inv_interp(0.0)

            # use math; speeds up because you don't need the full interpolated object
            for i_lz_p = 1:(G.N_lz-1)
                if (w_over_w_bar[i_lz_p] < 0) & (w_over_w_bar[i_lz_p+1] > 0)
                    m = (w_over_w_bar[i_lz_p+1] - w_over_w_bar[i_lz_p])/(G.grid_lz[i_lz_p+1] - G.grid_lz[i_lz_p])
                    b = w_over_w_bar[i_lz_p] - m * G.grid_lz[i_lz_p]
                    R.lz_d[i_k_p, i_b_p, i_lz] = -b/m
                    break
                end
            end
        end
    end
end

function compute_q!(P::Primitives, G::Grids, R::Results)
   
    q_next = zeros(G.N_k, G.N_b, G.N_lz) .+ 1/(1+P.r)

    for (i_k_p, k_p) = enumerate(G.grid_k), (i_b_p, b_p) = enumerate(G.grid_b), (i_lz, lz) = enumerate(G.grid_lz)
        if R.lz_d[i_k_p, i_b_p, i_lz] != G.min_lz
            recovery_return = 0.0
            nondefault_mass = 0.0
            for (i_lz_p, lz_p) = enumerate(G.grid_lz)
                if lz_p >= R.lz_d[i_k_p, i_b_p, i_lz]
                    nondefault_mass += G.Π_lz[i_lz, i_lz_p]
                else
                    recovery_value = (1-P.ξ)*(1-P.δ)*k_p + exp(lz_p) * k_p^P.α - T_C(exp(lz_p) * k_p^P.α - P.δ * k_p, P) - R.w_bar[i_lz_p]
                    if recovery_value < 0.0
                        println(R.w_bar[i_lz_p])
                        error("recovery value negative")
                    end
                recovery_return += recovery_value / b_p * G.Π_lz[i_lz, i_lz_p]
                
                end
            end
            # compute interest rate
            r_tilde_next = (1/(1-P.τ_i)) * ((1 + P.r * (1 - P.τ_i) - recovery_return)/nondefault_mass - 1)

            # convert to bond price
            q_next[i_k_p, i_b_p, i_lz] = max(min(1/(1+r_tilde_next), 1/(1+P.r)), 0.0)
        end
    end
    return q_next
end

function solve_model(;show_progress = true)
    P = Initialize_Primitives()
    G = Initialize_Grids(P)
    R = Initialize_Results(P, G)
    
    tolerence = 0.3
    err = 100.0
    i = 1
    max_iterations = 100

    λ = 0.95

    while true
        solve_vf!(P, G, R; show_progress = false)
        compute_default_states!(P, G, R)
        q_guess = compute_q!(P, G, R)
        err = maximum(abs.(R.q - q_guess))
        R.q = R.q .* λ .+ q_guess .* (1-λ);
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
        i+=1
    end
    return R
end