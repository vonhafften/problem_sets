####################################################################################################################
# ECON 899A Computational Economics
# Problem set 5

# Krusell and Smith (1998)
# Alex von Hafften
# October 13, 2021

# Based on code written by Phil Coyle

# This file contains the function for value function iteration.
####################################################################################################################

include("02_initialize.jl");

using Interpolations, Optim, Distributed

# Utility function
function u(c::Float64)
    if c > 0
       utility = log(c)
    else
        utility = -1/eps()
    end
    
    return utility
end

# Apply Bellman operator
function Bellman(results::Results, primitives::Primitives, grids::Grids, shocks::Shocks)
    @unpack β, α, δ, e_bar = primitives
    @unpack n_k, n_ε, n_K, n_z = grids
    @unpack k_min, k_max, K_min, K_max = grids
    @unpack k_grid, k_grid_srl, ε_grid, K_grid, K_grid_srl, z_grid = grids
    @unpack u_g, u_b, Π_ε = shocks

    @unpack value_function, policy_function, a0, a1, b0, b1 = results

    # Interpolate value function
    value_interp = interpolate(value_function, BSpline(Linear()))
    value_interp = Interpolations.scale(value_interp, k_grid_srl, 1:2, K_grid_srl, 1:2)

    # create updated value and policy functions
    v_next  = zeros(n_k, n_ε, n_K, n_z)
    pf_next = zeros(n_k, n_ε, n_K, n_z)

    # Iterate over aggregate states
    for (i_z, z_today) in enumerate(z_grid)

        # Set labor for today based on unemployment rate.
        if i_z == 1
            L_today = e_bar * (1 - u_g)
        elseif i_z == 2
            L_today = e_bar * (1 - u_b)
        end

        # Iterate over capital levels today
        for (i_K, K_today) in enumerate(K_grid)

            # interest rates and wages today 
            w_today = (1 - α) * z_today * (K_today/L_today) ^ α
            r_today = α * z_today * (K_today/L_today) ^ (α - 1)

            # Forecast aggregate capital
            if i_z == 1
                K_tomorrow = exp(a0 + a1 * log(K_today))
            else # elseif i_z == 2
                K_tomorrow = exp(b0 + b1 * log(K_today))
            end
            K_tomorrow = max(min(K_max, K_tomorrow), K_min)

            # Iterate over employment status today
            for (i_ε, ε_today) in enumerate(ε_grid)
                row = i_ε + n_ε * (i_z-1)

                # Iterate over individual capital holdings today
                for (i_k, k_today) in enumerate(k_grid)
                    budget = r_today * k_today + w_today * ε_today + (1.0 - δ) * k_today

                    # Define function of continuation value 
                    v_tomorrow(kp) = Π_ε[row,1]*value_interp(kp,1,K_tomorrow,1) + Π_ε[row,2]*value_interp(kp,2,K_tomorrow,1) +
                                     Π_ε[row,3]*value_interp(kp,1,K_tomorrow,2) + Π_ε[row,4]*value_interp(kp,2,K_tomorrow,2)

                    # Define function of functional equation
                    val_func(kp) = u(budget - kp) +  β * v_tomorrow(kp)

                    # Define negative (built-in method finds minimum)
                    obj(kp) = -val_func(kp)

                    # Use in-built Brent's method
                    opt = optimize(obj, k_min, min(k_max, budget))

                    # Update
                    v_next[i_k, i_ε, i_K, i_z] = -opt.minimum
                    pf_next[i_k, i_ε, i_K, i_z] = opt.minimizer
                end
            end
        end
    end
    return v_next, pf_next
end

function Solve_Bellman(results::Results)
    @unpack max_iterations, tol_vfi = Simulation()

    primitives = Primitives()
    grids = Grids()
    shocks = Shocks()

    error, i = 100, 1 # convergence variables

    while error > tol_vfi # loops until convergence
        
        # Applies Bellman operator
        v_next, pf_next = Bellman(results, primitives, grids, shocks)
        error = maximum(abs.(v_next .- results.value_function)) # sup norm

        results.value_function = v_next # update
        results.policy_function = pf_next # update

        if i >= max_iterations
            error("Value function did not converge.")
            break
        end
        i += 1
    end

    println("Value function converged in ", i, " iterations.")
    
    return results
end