#############################################################################
# Krusell and Smith (1998)
# Alex von Hafften
# October 13, 2021

# ECON 899A Computational Economics
# Problem set 5

# Based on code written by Phil Coyle

# This file contains the functions to estimate the model.
# ./run.jl calls this file and runs the model.
#############################################################################

using Parameters, Interpolations, Optim

## Housekeeping
@with_kw struct Params
    cBET::Float64 = 0.99
    cALPHA::Float64 = 0.36
    cDEL::Float64 = 0.025
    cLAM::Float64 = 0.5

    N::Int64 = 5000
    T::Int64 = 11000
    burn::Int64 = 1000

    tol_vfi::Float64 = 1e-4
    tol_coef::Float64 = 1e-4
    tol_r2::Float64 = 1.0 - 1e-2
    maxit::Int64 = 10000
end

@with_kw struct Grids
    k_lb::Float64 = 0.001
    k_ub::Float64 = 20.0
    n_k::Int64 = 21
    k_grid::Array{Float64,1} = range(k_lb, stop = k_ub, length = n_k)

    K_lb::Float64 = 10.0
    K_ub::Float64 = 15.0
    n_K::Int64 = 11
    K_grid::Array{Float64,1} = range(K_lb, stop = K_ub, length = n_K)

    eps_l::Float64 = 0.0
    eps_h::Float64 = 0.3271
    n_eps::Int64 = 2
    eps_grid::Array{Float64,1} = [eps_h, eps_l]

    z_l::Float64 = 0.99
    z_h::Float64 = 1.01
    n_z::Int64 = 2
    z_grid::Array{Float64,1} = [z_h, z_l]
end

@with_kw  struct Shocks
    #parameters of transition matrix:
    d_ug::Float64 = 1.5 # Unemp Duration (Good Times)
    u_g::Float64 = 0.04 # Fraction Unemp (Good Times)
    d_g::Float64 = 8.0 # Duration (Good Times)
    u_b::Float64 = 0.1 # Fraction Unemp (Bad Times)
    d_b::Float64 = 8.0 # Duration (Bad Times)
    d_ub::Float64 = 2.5 # Unemp Duration (Bad Times)

    #transition probabilities for aggregate states
    pgg::Float64 = (d_g-1.0)/d_g
    pgb::Float64 = 1.0 - (d_b-1.0)/d_b
    pbg::Float64 = 1.0 - (d_g-1.0)/d_g
    pbb::Float64 = (d_b-1.0)/d_b

    #transition probabilities for aggregate states and staying unemployed
    pgg00::Float64 = (d_ug-1.0)/d_ug
    pbb00::Float64 = (d_ub-1.0)/d_ub
    pbg00::Float64 = 1.25*pbb00
    pgb00::Float64 = 0.75*pgg00

    #transition probabilities for aggregate states and becoming employed
    pgg01::Float64 = (u_g - u_g*pgg00)/(1.0-u_g)
    pbb01::Float64 = (u_b - u_b*pbb00)/(1.0-u_b)
    pbg01::Float64 = (u_b - u_g*pbg00)/(1.0-u_g)
    pgb01::Float64 = (u_g - u_b*pgb00)/(1.0-u_b)

    #transition probabilities for aggregate states and becoming unemployed
    pgg10::Float64 = 1.0 - (d_ug-1.0)/d_ug
    pbb10::Float64 = 1.0 - (d_ub-1.0)/d_ub
    pbg10::Float64 = 1.0 - 1.25*pbb00
    pgb10::Float64 = 1.0 - 0.75*pgg00

    #transition probabilities for aggregate states and staying employed
    pgg11::Float64 = 1.0 - (u_g - u_g*pgg00)/(1.0-u_g)
    pbb11::Float64 = 1.0 - (u_b - u_b*pbb00)/(1.0-u_b)
    pbg11::Float64 = 1.0 - (u_b - u_g*pbg00)/(1.0-u_g)
    pgb11::Float64 = 1.0 - (u_g - u_b*pgb00)/(1.0-u_b)

    # Markov Transition Matrix
    Mgg::Array{Float64,2} = [pgg11 pgg01
                            pgg10 pgg00]

    Mbg::Array{Float64,2} = [pgb11 pgb01
                            pgb10 pgb00]

    Mgb::Array{Float64,2} = [pbg11 pbg01
                            pbg10 pbg00]

    Mbb ::Array{Float64,2} = [pbb11 pbb01
                             pbb10 pbb00]

    markov::Array{Float64,2} = [pgg*Mgg pgb*Mgb
                                pbg*Mbg pbb*Mbb]
end

mutable struct Results
    pf_k::Array{Float64,4}
    pf_v::Array{Float64,4}

    a0::Float64
    a1::Float64
    b0::Float64
    b1::Float64

    R2::Array{Float64,1}
end

function draw_shocks(S::Shocks, N::Int64,T::Int64)
    @unpack pgg, pbb, Mgg, Mgb, Mbg, Mbb = S

    # Shock
    Random.seed!(12032020)
    dist = Uniform(0, 1)

    # Allocate space for shocks and initialize
    idio_state = zeros(N,T)
    agg_state = zeros(T)
    idio_state[:,1] .= 1
    agg_state[1] = 1

    for t = 2:T
        agg_shock = rand(dist)
        if agg_state[t-1] == 1 && agg_shock < pgg
            agg_state[t] = 1
        elseif agg_state[t-1] == 1 && agg_shock > pgg
            agg_state[t] = 2
        elseif agg_state[t-1] == 2 && agg_shock < pbb
            agg_state[t] = 2
        elseif agg_state[t-1] == 2 && agg_shock > pbb
            agg_state[t] = 1
        end

        for n = 1:N
            idio_shock = rand(dist)
            if agg_state[t-1] == 1 && agg_state[t] == 1
                p11 = Mgg[1,1]
                p00 = Mgg[2,2]

                if idio_state[n,t-1] == 1 && idio_shock < p11
                    idio_state[n,t] = 1
                elseif idio_state[n,t-1] == 1 && idio_shock > p11
                    idio_state[n,t] = 2
                elseif idio_state[n,t-1] == 2 && idio_shock < p00
                    idio_state[n,t] = 2
                elseif idio_state[n,t-1] == 2 && idio_shock > p00
                    idio_state[n,t] = 1
                end
            elseif agg_state[t-1] == 1 && agg_state[t] == 2
                p11 = Mgb[1,1]
                p00 = Mgb[2,2]

                if idio_state[n,t-1] == 1 && idio_shock < p11
                    idio_state[n,t] = 1
                elseif idio_state[n,t-1] == 1 && idio_shock > p11
                    idio_state[n,t] = 2
                elseif idio_state[n,t-1] == 2 && idio_shock < p00
                    idio_state[n,t] = 2
                elseif idio_state[n,t-1] == 2 && idio_shock > p00
                    idio_state[n,t] = 1
                end
            elseif agg_state[t-1] == 2 && agg_state[t] == 1
                p11 = Mbg[1,1]
                p00 = Mbg[2,2]

                if idio_state[n,t-1] == 1 && idio_shock < p11
                    idio_state[n,t] = 1
                elseif idio_state[n,t-1] == 1 && idio_shock > p11
                    idio_state[n,t] = 2
                elseif idio_state[n,t-1] == 2 && idio_shock < p00
                    idio_state[n,t] = 2
                elseif idio_state[n,t-1] == 2 && idio_shock > p00
                    idio_state[n,t] = 1
                end
            elseif agg_state[t-1] == 2 && agg_state[t] == 2
                p11 = Mbb[1,1]
                p00 = Mbb[2,2]

                if idio_state[n,t-1] == 1 && idio_shock < p11
                    idio_state[n,t] = 1
                elseif idio_state[n,t-1] == 1 && idio_shock > p11
                    idio_state[n,t] = 2
                elseif idio_state[n,t-1] == 2 && idio_shock < p00
                    idio_state[n,t] = 2
                elseif idio_state[n,t-1] == 2 && idio_shock > p00
                    idio_state[n,t] = 1
                end
            end
        end
    end

    return idio_state, agg_state
end

function Bellman(P::Params, G::Grids, S::Shocks, R::Results)
    @unpack cBET, cALPHA, cDEL = P
    @unpack n_k, k_grid, n_eps, eps_grid, eps_h, K_grid, n_K, n_z, z_grid = G
    @unpack u_g, u_b, markov = S
    @unpack pf_k, pf_v, a0, a1, b0, b1= R

    pf_k_up = zeros(n_k, n_eps, n_K, n_z)
    pf_v_up = zeros(n_k, n_eps, n_K, n_z)

    # In Julia, this is how we define an interpolated function.
    # Need to use the package "Interpolations".
    # (If you so desire, you can write your own interpolation function too!)
    k_interp = interpolate(k_grid, BSpline(Linear()))
    v_interp = interpolate(pf_v, BSpline(Linear()))

    for (i_z, z_today) in enumerate(z_grid)
        for (i_K, K_today) in enumerate(K_grid)
            if i_z == 1
                K_tomorrow = a0 + a1*log(K_today)
            elseif i_z == 2
                K_tomorrow = b0 + b1*log(K_today)
            end
            K_tomorrow = exp(K_tomorrow)

            # See that K_tomorrow likely does not fall on our K_grid...this is why we need to interpolate!
            i_Kp = get_index(K_tomorrow, K_grid)

            for (i_eps, eps_today) in enumerate(eps_grid)
                row = i_eps + n_eps*(i_z-1)

                for (i_k, k_today) in enumerate(k_grid)
                    budget_today = r_today*k_today + w_today*eps_today + (1.0 - cDEL)*k_today

                    # We are defining the continuation value. Notice that we are interpolating over k and K.
                    v_tomorrow(i_kp) = markov[row,1]*v_interp(i_kp,1,i_Kp,1) + markov[row,2]*v_interp(i_kp,2,i_Kp,1) +
                                        markov[row,3]*v_interp(i_kp,1,i_Kp,2) + markov[row,4]*v_interp(i_kp,2,i_Kp,2)


                    # We are now going to solve the HH's problem (solve for k).
                    # We are defining a function val_func as a function of the agent's capital choice.
                    val_func(i_kp) = log(budget_today - k_interp(i_kp)) +  cBET*v_tomorrow(i_kp)

                    # Need to make our "maximization" problem a "minimization" problem.
                    obj(i_kp) = -val_func(i_kp)
                    lower = 1.0
                    upper = get_index(budget_today, k_grid)

                    # Then, we are going to maximize the value function using an optimization routine.
                    # Note: Need to call in optimize to use this package.
                    opt = optimize(obj, lower, upper)

                    k_tomorrow = k_interp(opt.minimizer[1])
                    v_today = -opt.minimum

                    # Update PFs
                    pf_k_up[i_k, i_eps, i_K, i_z] = k_tomorrow
                    pf_v_up[i_k, i_eps, i_K, i_z] = v_today
                end
            end
        end
    end

    return pf_k_up, pf_v_up
end


@unpack T, N = Params()
@unpack n_k, n_K, n_eps, n_z, k_grid = Grids()

results = Results(zeros(n_k, n_eps, n_K, n_z), zeros(n_k, n_eps, n_K, n_z), 0.095, 0.999, 0.085, 0.999, zeros(n_k))


Bellman(Params(), Grids(), Shocks(), results)

function get_index(val::Float64, grid::Array{Float64,1})
    n = length(grid)
    index = 0
    if val <= grid[1]
        index = 1
    elseif val >= grid[n]
        index = n
    else
        index_upper = findfirst(x->x>val, grid)
        index_lower = index_upper - 1
        val_upper, val_lower = grid[index_upper], grid[index_lower]

        index = index_lower + (val - val_lower) / (val_upper - val_lower)
    end
    return index
end
