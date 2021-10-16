####################################################################################################################
# Krusell and Smith (1998)
# Alex von Hafften
# October 13, 2021

# ECON 899A Computational Economics
# Problem set 5

# This file contains runs the model.
####################################################################################################################

include("03_value_function.jl");

function simulate_capital_path(results::Results)
    @unpack k_grid_srl, K_grid_srl, n_ε, n_z = Grids()
    @unpack N, T, k_ss = Simulation()
    
    @unpack Z, E, policy_function, K = results

    # interpolate policy function
    policy_interp = interpolate(policy_function, BSpline(Linear()))
    policy_interp = Interpolations.scale(policy_interp, k_grid_srl, 1:n_ε, K_grid_srl, 1:n_z)

    # initialize matrix of capital over agents and time
    k_holdings = zeros(N, T)

    # first period starts at complete markets steady state
    k_holdings[:, 1] .= k_ss
    K[1] = k_ss

    for t = 2:T # iterate over time
        for i = 1:N # iterate over agents
            # get next periods capital holdings from policy funciton
            k_holdings[i, t] = policy_interp(k_holdings[i, t-1], E[i, t-1], K[t-1], Z[t-1])
        end

        # compute aggregate capital 
        K[t] = sum(k_holdings[:, t])/N  
    end

    return K
end

function estimate_regression(results::Results)
    @unpack burn_in, T = Simulation()
    @unpack Z, K = results

    # remove burn-in period
    Kp = K[burn_in:T]
    K  = K[burn_in-1:T-1]
    Z  = Z[burn_in-1:T-1]

    # boom periods
    y_a = log.(Kp[Z .== 1])
    x_a = [ones(length(y_a)) log.(K[Z .== 1])]

    beta_a = inv(x_a' * x_a) * x_a' * y_a

    ssr_a = sum((y_a - x_a * beta_a).^2)
    sst_a = sum((y_a .- mean(y_a)).^2)

    # recession periods
    y_b = log.(Kp[Z .== 2])
    x_b = [ones(length(y_b)) log.(K[Z .== 2])]

    beta_b = inv(x_b' * x_b) * x_b' * y_b

    ssr_b = sum((y_b - x_b * beta_b).^2)
    sst_b = sum((y_b .- mean(y_b)).^2)

    R2 = 1 - (ssr_a + ssr_b)/(sst_a + sst_b)

    return vcat(beta_a, beta_b, R2)
end

function Solve_model()
    @unpack max_iterations, tol_coef, tol_R2, λ = Simulation()

    results = Initialize()

    err_coef, err_R2, i = 100, 100, 1

    println("********************************************************")

    while err_coef > tol_coef || results.R2 < tol_R2

        println("Iteration #", i)

        @unpack a0, a1, b0, b1, R2 = results

        # VFI, simulates capital, and estimates regression
        results = Solve_Bellman(results)
        results.K = simulate_capital_path(results)
        a0_hat, a1_hat, b0_hat, b1_hat, R2_hat = estimate_regression(results)

        err_coef = abs(a0 - a0_hat) + abs(a1 - a1_hat) + abs(b0 - b0_hat) + abs(b1 - b1_hat)

        results.a0 = λ * a0_hat + (1 - λ) * a0
        results.a1 = λ * a1_hat + (1 - λ) * a1
        results.b0 = λ * b0_hat + (1 - λ) * b0
        results.b1 = λ * b1_hat + (1 - λ) * b1

        results.R2 = R2_hat

        println("Coefficient error is ", err_coef)
        println("R-squared is ", R2_hat)
        println("********************************************************")

        if i >= max_iterations
            println("Coefficients did not converge.")
            break
        end

        i += 1
    end

    println("Coefficients converged")

    results

end