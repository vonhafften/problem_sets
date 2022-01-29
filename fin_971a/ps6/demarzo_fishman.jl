# Problem set 6 of FIN 971
# Professor: Dean Corbae
# DeMarzo Fishman (2007, RFS)

# Alex von Hafften
# December 14, 2021

using Parameters, Plots, JuMP, GLPK

cd("/Users/alexandervonhafften/Documents/UW Madison/problem_sets/fin_971/ps6/")

# struct to hold model primitives (global constants)
@with_kw struct Primitives
    β::Float64               = exp(-0.0953)                               # principal discount factor
    δ::Float64               = exp(-0.0998)                               # agent discount factor
    L::Float64               = 75.0                                       # liquidation value for principal
    R::Float64               = 0.0                                        # recovery value for agent
    λ::Float64               = 1.0                                        # moral hazard parameter
    Y_L                      = 0.0                                        # low cash-flow
    Y_H                      = 20.0                                       # high cash-flow
    π_L::Float64             = 0.5                                        # probability of low cash flow
    π_H::Float64             = 0.5                                        # probability of high cash flow
    μ::Float64               = π_H * Y_H + π_L * Y_L                      # mean cash flow
    a_min::Float64           = max(π_H * λ * (Y_H - Y_L) + R, R/δ) + 1e-1 # minimum promised utility level (eq. 13 in handout on DF)
    # a_max::Float64           = μ/(1 - β)                                # maximum promised utility level (eq. 14 in handout on DF)
    a_max::Float64           = 60.0                                       # maximum promised utility level (eq. 14 in handout on DF)
    a_n::Int64               = 101                                        # number of promised utility level states
    a_grid                   = range(a_min, a_max; length = a_n)          # promised utility level grid
end

# struct to hold model outputs
mutable struct Results
    b::Array{Float64,1}      # value  - expected utility of investor
    p_L::Array{Float64,1}    # policy - probability of termination after low report
    p_H::Array{Float64,1}    # policy - probability of termination after high report
    d_L::Array{Float64,1}    # policy - equity injection after low report
    d_H::Array{Float64,1}    # policy - equity injection after high report
    a_p_L::Array{Float64,1}  # policy - promized future utility to agent after low report
    a_p_H::Array{Float64,1}  # policy - promized future utility to agent after high report
end

# first-best for initial guess; see section 2.3 of handout on DF
function b_fb()
    @unpack a_grid, μ, β, δ, R = Primitives()
    
    K_0 = (μ - (β - δ)/δ)/(1 - β)
    μ + β * K_0 - (β - δ) * R / δ .- a_grid 
end

# initialize results structure
function Initialize() 
    @unpack a_n = Primitives()

    Results(b_fb(), zeros(a_n), zeros(a_n), zeros(a_n), zeros(a_n), zeros(a_n), zeros(a_n))
end

# Bellman operator
function Bellman(results::Results, primitives::Primitives)
    @unpack π_H, β, δ, L, R, λ, Y_L, Y_H, π_L, π_H, μ, a_min, a_max, a_n, a_grid = primitives # unpack primitive structure
    
    # initialize structure to store updated results in.
    results_next = Initialize() 

    i_a_p_L_start = 1
    i_a_p_H_start = 1
    previous_max = 100

    # loop over state space
    for (i, a) = enumerate(a_grid)

        candidate_max = -1/eps()

        # grid search over the future promised utility levels
        for i_a_p_L = i_a_p_L_start:a_n

            # we know that a high report will mean higher promised utility so start no lower than the low promised level
            for i_a_p_H = max(i_a_p_L, i_a_p_H_start):a_n

                # get promised levels
                a_p_L = a_grid[i_a_p_L]
                a_p_H = a_grid[i_a_p_H]
            
                # uses JuMP for constrained optimization over p_L, p_H, d_L, and d_H
                model = Model(GLPK.Optimizer)

                # add variables
                @variable(model, 0 <= p_L <= 1)
                @variable(model, 0 <= p_H <= 1)
                @variable(model, 0 <= d_L)
                @variable(model, 0 <= d_H)
        
                # add objective
                @objective(model, Max, π_H*(Y_H - d_H + p_H*L) + π_L*(Y_L - d_L + p_L *L) + β*(π_H*(1-p_H)*results.b[i_a_p_H] + π_L*(1-p_L)*results.b[i_a_p_L]))

                # add constraints
                @constraint(model, c1, a == π_H*(d_H + p_H*R) + π_L * (d_L + p_L*R) + δ * (π_H*(1-p_H)*a_p_H + π_L *(1-p_L)*a_p_L))
                @constraint(model, c2, d_H + (1-p_H)*δ*a_p_H + p_H*R == λ*(Y_H-Y_L) + d_L + (1-p_L)*δ*a_p_L + p_L*R )
                @constraint(model, c3, p_H*R + (1-p_H)*δ*a_p_H >= R)
                @constraint(model, c4, p_L*R + (1-p_L)*δ*a_p_L >= R)

                # optimization step
                optimize!(model)

                # if optimized level higher than current update results_next structure
                if (objective_value(model) > candidate_max)
                    
                    results_next.b[i] = objective_value(model)

                    results_next.p_L[i] = value(p_L)
                    results_next.p_H[i] = value(p_H)
                    results_next.d_L[i] = value(d_L)
                    results_next.d_H[i] = value(d_H)
                    results_next.a_p_L[i] = a_p_L
                    results_next.a_p_H[i] = a_p_H

                    i_a_p_L_start = i_a_p_L
                    i_a_p_H_start = i_a_p_H
                    candidate_max = objective_value(model)
                else # if we get to a point that the max is not begin updated then we can break
                    break
                end
            end
        end        
    end
    results_next
end

# sup norm over value and policy functions
function sup_norm(results_1::Results, results_2::Results)
    maximum(vcat(abs.(results_1.b .- results_2.b), 
                 abs.(results_1.p_H .- results_2.p_H),
                 abs.(results_1.p_L .- results_2.p_L), 
                 abs.(results_1.d_L .- results_2.d_L), 
                 abs.(results_1.d_H .- results_2.d_H), 
                 abs.(results_1.a_p_L .- results_2.a_p_L), 
                 abs.(results_1.a_p_H .- results_2.a_p_H)))
end

# function to solve the model
function Solve_model()

    primitives = Primitives()
    results = Initialize() 

    # loop until convergence
    error, i, max_iter = 100, 0, 100
    while (error > 0.001) & (i < max_iter)
        i += 1
        results_next = Bellman(results, primitives)         # next guess of b function
        error = sup_norm(results, results_next) # check for convergence
        results = results_next                             # update
        println("Current error: ", error)
    end

    println("Value function converged in ", i, " iterations")
    results
end

primitives = Primitives();
@elapsed results = Solve_model()

# Value function plot
plot(primitives.a_grid, results.b, label = "Assymetric Info", title = "Value Function");
plot!(primitives.a_grid, b_fb(), label = "First Best");
plot!(size=(400,400));
savefig("value_function.png")

# Policy function plot
p1 = plot(primitives.a_grid, [results.a_p_H primitives.a_grid], ylim = (10, 60), title = "a'_H");
p2 = plot(primitives.a_grid, [results.a_p_L primitives.a_grid], ylim = (10, 60), title = "a'_L");
p3 = plot(primitives.a_grid, results.p_H, ylim = (0, 1), title = "p_H");
p4 = plot(primitives.a_grid, results.p_L, ylim = (0, 1), title = "p_L");
p5 = plot(primitives.a_grid, results.d_H, ylim = (0, 20), title = "d_H");
p6 = plot(primitives.a_grid, results.d_L, ylim = (0, 20), title = "d_L");

plot(p1, p2, p3, p4, p5, p6, layout =  (3, 2), legend = false);
plot!(size=(400,600));
savefig("policy_functions.png")

