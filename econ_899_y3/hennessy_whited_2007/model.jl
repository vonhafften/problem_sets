# Code to solve Hennessy and Whited 2007

# Alex von Hafften
# Sept. 27, 2022

using Parameters, QuantEcon, Interpolations, Optim

include("structures.jl")
include("helper_functions.jl")

# applies bellman operator
function Apply_Bellman!(P::Primitives, G::Grids, R::Results, grid_search::Bool)

    # interpolate interest rates and value
    r_tilde_interp = LinearInterpolation((G.grid_k, G.grid_b, G.grid_lz_c), R.r_tilde)
    vf_interp = LinearInterpolation((G.grid_w, G.grid_lz_c), R.vf)

    # initialize next guess for value function
    vf_next = zeros(G.N_w, G.N_lz_c)

    # iterate over state variables
    for (i_w, w) = enumerate(G.grid_w), (i_lz, lz) = enumerate(G.grid_lz_c)
        
        # function to compute payoff
        function payoff(k_p::Float64, b_p::Float64)
            # determine whether cash distribution or equity issuance
            Φ = (w + b_p - k_p > 0)

            # flow payoff
            result = Φ * (w + b_p - k_p - compute_T_d(w + b_p - k_p, P))
            result -= (1-Φ) * (k_p - w - b_p + compute_Λ(k_p - w - b_p, P))

            # determine interest rate
            r_tilde = r_tilde_interp(k_p, b_p, lz)

            # add continuation value
            for (i_lz_p, lz_p) = enumerate(G.grid_lz_c)
                w_p = compute_nw(k_p, b_p, exp(lz), exp(lz_p), r_tilde, P)
                result += 1/(1+P.r*(1 - P.τ_i)) * G.Π_lz_c[i_lz, i_lz_p] * max(vf_interp(w_p, lz_p), 0.0)
            end

            return result
        end

        # solve with grid search
        if (grid_search)
            candidate_max = -1/eps()
            for (i_k_p, k_p) = enumerate(G.grid_k), (i_b_p, b_p) = enumerate(G.grid_b)
                candidate = payoff(k_p, b_p)
                if (candidate > candidate_max)
                    candidate_max = candidate
                    vf_next[i_w, i_lz] = candidate
                    R.pf_b[i_w, i_lz] = b_p
                    R.pf_k[i_w, i_lz] = k_p
                end
            end
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


    end
    return vf_next
end


# solves for policy and value function with bond prices given
function Solve_vf!(P::Primitives, G::Grids, R::Results)
    err = 100
    tolerence = 0.001
    i = 1
    max_iterations = 10
    grid_search = true
    while true
        println(err)
        vf_next = Apply_Bellman!(P, G, R, grid_search)
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

