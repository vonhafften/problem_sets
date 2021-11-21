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
function Φ(x::Float64)
    return cdf(Normal(0, 1), x)
end

# inverse standard normal cdf
function Φ_inverse(p::Float64)
    return quantile(Normal(0, 1), p)
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
        result = Φ((α_0 + x'*β + z'*γ)/σ_0)
    elseif t == 2.0
        f_2(ε_0) = Φ(α_1 + x'*β + z'*γ + ρ*ε_0) * ϕ(ε_0/σ_0) / σ_0
        result = integrate_1d(f_2, α_0 + x'*β + z'*γ, KPU_1d)
    elseif t == 3.0
        f_3(ε_0, ε_1) = Φ(-α_2 - x'*β - z'*γ - ρ*ε_1) * ϕ(ε_1 - ρ*ε_0) * ϕ(ε_0/σ_0) / σ_0
        result = integrate_2d(f_3, α_0 + x'*β + z'*γ, α_1 + x'*β + z'*γ, KPU_2d)
    elseif t == 4.0
        f_4(ε_0, ε_1) = Φ( α_2 + x'*β + z'*γ - ρ*ε_1) * ϕ(ε_1 - ρ*ε_0) * ϕ(ε_0/σ_0) / σ_0
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
    t::Float64, x::Array{Float64, 1}, z::Array{Float64, 1}, u_0::Array{Float64, 1}, u_1::Array{Float64, 1}, u_2::Array{Float64, 1})

    n_trials = length(u_0)
    σ_0 = 1/(1 - ρ)^2

    truncation_0 = Φ((-α_0 - x'*β - z'*γ)/σ_0) # evaluates truncation point for first shock probability

    if t == 1.0

        return 1 - truncation_0
    
    else # if t = 2.0 or 3.0 or 4.0
    
        pr_0 = u_0 * truncation_0 # scales uniform rv between zero and the truncation point.
        η_0 = Φ_inverse.(pr_0)
        ε_0 = η_0 .* σ_0

        truncation_1 = Φ.(-α_1 .- x'*β .- z'*γ .- ρ.*ε_0) # initializes simulation-specific truncation points

        if t == 2.0

            return sum(truncation_0 .* (1 .- truncation_1))/n_trials
    
        else # if t = 3.0 or 4.0
    
            pr_1 = u_1 .* truncation_1 # scales uniform rv between zero and the truncation point.
            η_1 = Φ_inverse.(pr_1) # initializes first shocks
            ε_1 = ρ .* ε_0 .+ η_1 

            truncation_2 = Φ.(-α_2 .- x'*β .- z'*γ .- ρ.*ε_1)
    
            if t == 3.0
                
                return sum(truncation_0 .* truncation_1 .* (1 .- truncation_2))/n_trials
    
            else # t = 4.0
    
                return sum(truncation_0 .* truncation_1 .* truncation_2)/n_trials

            end
        end
    end
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
        count = sum( α_0 + x'*β + z'*γ .+ ε_0 .> 0)
    elseif t == 2.0
        count = sum((α_0 + x'*β + z'*γ .+ ε_0 .< 0) .* (α_1 + x'*β + z'*γ .+ ε_1 .> 0))
    elseif t == 3.0
        count = sum((α_0 + x'*β + z'*γ .+ ε_0 .< 0) .* (α_1 + x'*β + z'*γ .+ ε_1 .< 0) .* (α_2 + x'*β + z'*γ .+ ε_2 .> 0))
    elseif t == 4.0
        count = sum((α_0 + x'*β + z'*γ .+ ε_0 .< 0) .* (α_1 + x'*β + z'*γ .+ ε_1 .< 0) .* (α_2 + x'*β + z'*γ .+ ε_2 .< 0))
    else
        error("Invalid value of t.")
    end

    # returns the frequency of the accepted simulations
    return count/length(ε_0)
end
