# Computational Economics
# Professor JF Houde
# Problem set 2
# Alex von Hafften 
# November 22, 2021

# Contains all function to go to worker processes

using Optim, ProgressMeter, DataFrames, Distributions

############################################################
# one-dimensional quadrature integration
############################################################

function integrate_1d(f, upper_bound, KPU_1d)

    # define functions to translate the (0, 1) interval into appropriate interval
    points = log.(KPU_1d[:, :x1]) .+ upper_bound
    jacobian = 1 ./KPU_1d[:, :x1]

    # sum over grid points
    return sum(KPU_1d[:, :weight] .* f.(points) .* jacobian)

end

############################################################
# two-dimensional quadrature integration
############################################################

function integrate_2d(f, upper_bound_0::Float64, upper_bound_1::Float64, KPU_2d)

    points_0 = log.(KPU_2d[:, :x1]) .+ upper_bound_0
    jacobian_0 = 1 ./ KPU_2d[:, :x1]

    points_1 = log.(KPU_2d[:, :x2]) .+ upper_bound_1
    jacobian_1 = 1 ./KPU_2d[:, :x2]

    return sum(KPU_2d[:, :weight] .* f.(points_0, points_1) .* jacobian_0 .* jacobian_1)

end

############################################################
# standard normal distribution functions
############################################################

# standard normal pdf
function ϕ(x::Float64)
    1/sqrt(2 * π) * exp((-1/2)*x^2)
end

# standard normal cdf
function Φ(x::Float64, KPU_1d; method = "built_in")
    if method == "built_in"

        if x > 4.0 # for really large values, the quadrature integration was resulting in probability above 1. By trial and error just return 1 at 4.0 or higher works.
            return 1.0
        else
            return integrate_1d(ϕ, x, KPU_1d)
        end
    elseif method == "quadrature"
        return cdf(Normal(0, 1), x)
    else
        error("Issue in Φ function")
    end
end

# inverse standard normal cdf
function Φ_inverse(p::Float64, KPU_1d; method = "built_in")
    if p < 0
        error("Invalid input to Φ_inverse function.")
    elseif p > 1
        error("Invalid input to Φ_inverse function.")
    elseif method == "custom"
        f(q) = abs(Φ(q[1], KPU_1d) - p)
        opt = optimize(f, [0.0,], GradientDescent())
        return opt.minimizer[1]
    elseif method == "built_in"
        return quantile(Normal(0, 1), p)
    else
        error("Issue in Φ_inverse function.")
    end
end

##############################################################################
# functions to calculate likelihood using quadrature integration methods
##############################################################################

# using quadrature integration
function likelihood_quadrature(γ::Array{Float64, 1}, β::Array{Float64, 1}, ρ::Float64, α_0::Float64, α_1::Float64, α_2::Float64, 
    t::Float64, x::Array{Float64, 1}, z::Array{Float64, 1}, KPU_1d, KPU_2d)
    
    result = 0.0
    σ_0 = 1/(1 - ρ)^2

    if t == 1.0
        result = Φ((-α_0 - x'*β - z'*γ)/σ_0, KPU_1d)
    elseif t == 2.0
        f_2(ε_0) = Φ(-α_1 - x'*β - z'*γ - ρ*ε_0, KPU_1d) * ϕ(ε_0/σ_0) / σ_0
        result = integrate_1d(f_2, α_0 + x'*β + z'*γ, KPU_1d)
    elseif t == 3.0
        f_3(ε_0, ε_1) = Φ(-α_2 - x'*β - z'*γ - ρ*ε_1, KPU_1d) * ϕ(ε_1 - ρ*ε_0) * ϕ(ε_0/σ_0) / σ_0
        result = integrate_2d(f_3, α_0 + x'*β + z'*γ, α_1 + x'*β + z'*γ, KPU_2d)
    elseif t == 4.0
        f_4(ε_0, ε_1) = Φ( α_2 + x'*β + z'*γ - ρ*ε_1, KPU_1d) * ϕ(ε_1 - ρ*ε_0) * ϕ(ε_0/σ_0) / σ_0
        result = integrate_2d(f_4, α_0 + x'*β + z'*γ, α_1 + x'*β + z'*γ, KPU_2d)
    else
        error("Invalid value of t.")
    end

    return result
end

function jacobian_quadrature(γ::Array{Float64, 1}, β::Array{Float64, 1}, ρ::Float64, α_0::Float64, α_1::Float64, α_2::Float64, 
    t::Float64, x::Array{Float64, 1}, z::Array{Float64, 1}, KPU_1d, KPU_2d)

    result = zeros(length(γ) + length(β) + 4)
    σ_0 = 1/(1 - ρ)^2

    if t == 1.0
        result = ϕ((-α_0 - x'*β - z'*γ)/σ_0) .* vcat(-z./σ_0, -x./ σ_0, -2*(1-ρ), -1, 0, 0)
    elseif t == 2.0
        ##################################
        
    elseif t == 3.0
        ##################################
        

    elseif t == 4.0
        ##################################
        

    else
        ##################################
        
    end

    return result

end


##############################################################################
# functions to calculate likelihood using GHK simulations
##############################################################################

# using ghk simulation
function likelihood_ghk(γ::Array{Float64, 1}, β::Array{Float64, 1}, ρ::Float64, α_0::Float64, α_1::Float64, α_2::Float64, 
    t::Float64, x::Array{Float64, 1}, z::Array{Float64, 1}, KPU_1d, KPU_2d, u_0::Array{Float64, 1}, u_1::Array{Float64, 1}, u_2::Array{Float64, 1})

    n_trials = length(u_0)
    
    result = 0.0
    σ_0 = 1/(1 - ρ)^2

    if t == 1.0

        truncation_0 = Φ((-α_0 - x'*β - z'*γ)/σ_0, KPU_1d) # evaluates truncation point for first shock probability
        pr_0 = u_0 * truncation_0 # scales uniform rv between zero and the truncation point
        result = sum(pr_0)/n_trials # result is average probability
        result_1 = sum(pr_0)/n_trials # result is average probability
    
    else # if t = 2.0 or 3.0 or 4.0
    
        truncation_0 = Φ((α_0 + x'*β + z'*γ)/σ_0, KPU_1d) # evaluates truncation for first shock (i.e. complement of trunctation point for t=1.0)
        pr_0 = u_0 * truncation_0 # scales uniform rv between zero and the truncation point.
        ε_0 = zeros(n_trials) # initializes first shocks
    
        for i = 1:n_trials # computes first shocks based on scaled probability
            ε_0[i] = Φ_inverse(pr_0[i], KPU_1d) * σ_0
        end
    
        if t == 2.0
    
            truncation_1 = zeros(n_trials) # initializes simulation-specific truncation points
            for i = 1:n_trials # calculates simulation-specific truncation points for second shock based on realization of first shock
                truncation_1[i] = Φ(-α_1 - x'*β - z'*γ - ρ * ε_0[i], KPU_1d)
            end
            pr_1 = u_1 .* truncation_1 # scales uniform rv random between zero and the truncation point.
            # result = sum(pr_0.*pr_1)/n_trials # result is average probability
            result_2 = sum(pr_1)/n_trials # result is average probability
    
        else # if t = 3.0 or 4.0
    
            truncation_1 = zeros(n_trials)
            for i = 1:n_trials
                truncation_1[i] = Φ(α_1 + x'*β + z'*γ + ρ * ε_0[i], KPU_1d) # complement of truncation point for t = 2.0
            end
            pr_1 = u_1 .* truncation_1 # scales uniform rv between zero and the truncation point.
            ε_1 = zeros(n_trials) # initializes first shocks
            for i = 1:n_trials # computes first shocks based on scaled probability and first shock realization
                ε_1[i] = Φ_inverse(pr_1[i], KPU_1d) + ρ*ε_0[i] # fills 
            end
    
            if t == 3.0
    
                truncation_2 = zeros(n_trials)
                for i = 1:n_trials # calculates simulation-specific truncation points for second shock based on realization of first shock
                    truncation_2[i] = Φ(-α_1 - x'*β - z'*γ - ρ * ε_1[i], KPU_1d)
                end
                pr_2 = u_2 .* truncation_2 # scales uniform rv random between zero and the truncation point.
                # result = sum(pr_0 .* pr_1 .* pr_2)/n_trials # result is average probability
                result_3 = sum(pr_2)/n_trials # result is average probability
    
            else # t = 4.0
    
                truncation_2 = zeros(n_trials)
                for i = 1:n_trials
                    truncation_2[i] = Φ(α_1 + x'*β + z'*γ + ρ * ε_1[i], KPU_1d) # complement of truncation point for t = 2.0
                end
                pr_2 = u_2 .* truncation_2 # scales uniform rv between zero and the truncation point.
                # result = sum(pr_0 .* pr_1 .* pr_2)/n_trials # result is average probability
                result_4 = sum(pr_2)/n_trials # result is average probability
    
            end
        end
    end
        
    return result
end

##############################################################################
# functions to calculate likelihood using accept-reject simulation
##############################################################################

function likelihood_accept_reject(γ::Array{Float64, 1}, β::Array{Float64, 1}, ρ::Float64, α_0::Float64, α_1::Float64, α_2::Float64, 
    t::Float64, x::Array{Float64, 1}, z::Array{Float64, 1}, ε_0::Array{Float64, 1}, ε_1::Array{Float64, 1}, ε_2::Array{Float64, 1})
    
    # initialize count variable
    count = 0

    # based on the value of t counts the number of accepted simulations
    if t == 1.0
        count = sum(ε_0 .< - α_0 - x'*β - z'*γ)
    elseif t == 2.0
        count = sum((ε_0 .< α_0 + x'*β + z'*γ) .*(ε_1 .< -α_1 .- x'*β .- z'*γ))
    elseif t == 3.0
        count = sum((ε_0 .< α_0 + x'*β + z'*γ) .*(ε_1 .< α_1 .+ x'*β .+ z'*γ) .* (ε_2 .< - α_2 .- x'*β .- z'*γ))
    elseif t == 4.0
        count = sum((ε_0 .< α_0 + x'*β + z'*γ) .*(ε_1 .< α_1 .+ x'*β .+ z'*γ) .* (ε_2 .< α_2 .+ x'*β .+ z'*γ))
    else
        error("Invalid value of t.")
    end

    # returns the frequency of the accepted simulations
    return count/length(ε_0)
end
