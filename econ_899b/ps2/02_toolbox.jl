# Computational Economics
# Professor JF Houde
# Problem set 2
# Alex von Hafften 
# November 22, 2021

# Contains all function to go to worker processes

using Optim, CSV, DataFrames

############################################################
# halton sequence generator
############################################################

# generates halton sequence of length n using base (base must be prime)
function halton(base::Int64, n::Int64)
    
    m, d = 0, 1
    burn_in = 500
    result = zeros(burn_in + n)

    # Adapted from code here: https://en.wikipedia.org/wiki/Halton_sequence
    for i = 1:(burn_in + n)
        x = d - m
        if x == 1
            m = 1
            d *= base
        else
            y = d / base
            while x <= y
                y /= base
            end
            m = (base + 1) * y - x
        end
        result[i] = m / d
    end

    result[(burn_in + 1):(burn_in + n)]

end

############################################################
# one-dimensional quadrature integration
############################################################

# quadrature nodes and weights
KPU_1d = DataFrame(CSV.File("PS2/KPU_d1_l20.csv"))

function integrate_1d(f, upper_bound)

    # define functions to translate the (0, 1) interval into appropriate interval
    points = log.(KPU_1d[:, :x1]) .+ upper_bound
    jacobian = 1 ./KPU_1d[:, :x1]

    # sum over grid points
    return sum(KPU_1d[:, :weight] .* f.(points) .* jacobian)

end

############################################################
# two-dimensional quadrature integration
############################################################

# quadrature nodes and weights
KPU_2d = DataFrame(CSV.File("PS2/KPU_d2_l20.csv"))

function integrate_2d(f, upper_bound_0, upper_bound_1)

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
    integrate_1d(ϕ, x)
end

# inverse standard normal cdf
function Φ_inverse(p::Float64)
    if p < 0
        error("Invalid input to Φ_inverse function.")
    elseif p > 1
        error("Invalid input to Φ_inverse function.")
    else
        f(q) = abs(Φ(q[1]) - p)
        opt = optimize(f, [0.0,], GradientDescent())
        return opt.minimizer[1]
    end
end

##############################################################################
# functions to calculate observation-level likelihood using different methods
##############################################################################

# using quadrature integration
function likelihood_quadrature(γ::Array{Float64, 1}, β::Array{Float64, 1}, ρ::Float64, α_0::Float64, α_1::Float64, α_2::Float64, 
    t::Float64, x::Array{Float64, 1}, z::Array{Float64, 1})
    
    result = 0.0
    σ_0 = 1/(1 - ρ)^2

    if t == 1.0
        result = Φ((-α_0 - x'*β - z'*γ)/σ_0)

    elseif t == 2.0

        f_2(ε_0) = Φ(-α_1 - x'*β - z'*γ - ρ*ε_0) * ϕ(ε_0/σ_0) / σ_0
        result   = integrate_1d(f_2, α_0 + x'*β + z'*γ)
    
    elseif t == 3.0

        f_3(ε_0, ε_1) = Φ(-α_2 - x'*β - z'*γ - ρ*ε_1) * ϕ(ε_1 - ρ*ε_0) * ϕ(ε_0/σ_0) / σ_0
        result        = integrate_2d(f_3, α_0 + x'*β + z'*γ, α_1 + x'*β + z'*γ)
    
    elseif t == 4.0

        f_4(ε_0, ε_1) = Φ( α_2 + x'*β + z'*γ - ρ*ε_1) * ϕ(ε_1 - ρ*ε_0) * ϕ(ε_0/σ_0) / σ_0
        result        = integrate_2d(f_4, -α_0 - x'*β - z'*γ, -α_1 - x'*β - z'*γ)

    else
        error("Invalid value of t.")
    end
    return result
end


# using accept-reject simulation
function likelihood_accept_reject(γ::Array{Float64, 1}, β::Array{Float64, 1}, ρ::Float64, α_0::Float64, α_1::Float64, α_2::Float64, 
    t::Float64, x::Array{Float64, 1}, z::Array{Float64, 1})
    
    n_trials = 100
    σ_0 = 1/(1 - ρ)^2

    # pulls independent shocks
    ε_0 = Φ_inverse.(halton(5,  n_trials)) * σ_0
    η_1 = Φ_inverse.(halton(7,  n_trials))
    η_2 = Φ_inverse.(halton(11, n_trials))

    # calculates error terms
    ε_1 = ρ .* ε_0 .+ η_1
    ε_2 = ρ .* ε_1 .+ η_2
    
    # initialize count variable
    count = 0

    # based on the value of t counts the number of accepted simulations
    if t == 1.0
        count = sum(ε_0 .< - α_0 - x'*β - z'*γ)
    elseif t == 2.0
        count = sum((ε_0 .< α_0 + x'*β + z'*γ) .*(η_1 .< -α_1 .- x'*β .- z'*γ .- ρ .* ε_0))
    elseif t == 3.0
        count = sum((ε_0 .< α_0 + x'*β + z'*γ) .*(ε_1 .< α_1 .+ x'*β .+ z'*γ) .* (ε_2 .< - α_2 .- x'*β .- z'*γ))
    elseif t == 4.0
        count = sum((ε_0 .< α_0 + x'*β + z'*γ) .*(ε_1 .< α_1 .+ x'*β .+ z'*γ) .* (ε_2 .< α_2 .+ x'*β .+ z'*γ))
    else
        error("Invalid value of t.")
    end

    # returns the frequency of the accepted simulations
    return count/n_trials
end

# compute likelihood for the matrix
function likelihood(γ::Array{Float64, 1}, β::Array{Float64, 1}, ρ::Float64, α_0::Float64, α_1::Float64, α_2::Float64,
    t::Array{Float64, 2}, x::Array{Float64, 2}, z::Array{Float64, 2})
    
    N = size(x)[1]
    
    # uses pmap; on a test with quadrature, this took about 3 minutes for the dataset
    # result = pmap(i-> likelihood(γ, β, ρ, α_0, α_1, α_2, t[i], x[i,:], z[i,:]), i for i = 1:N)
    
    # uses distributed for loop; on a test with quadrature, this took about 2 minutes for the dataset
    result = SharedArray{Float64}(N)
    @sync @distributed for i = 1:N
        result[i] = likelihood_ar(γ, β, ρ, α_0, α_1, α_2, t[i], x[i,:], z[i,:])
    end

    return result
end

# compute likelihood for the matrix
function likelihood(γ::Array{Float64, 1}, β::Array{Float64, 1}, ρ::Float64, α_0::Float64, α_1::Float64, α_2::Float64,
    t::Array{Float64, 2}, x::Array{Float64, 2}, z::Array{Float64, 2}; method = "quadrature")
    
    N = size(x)[1]
    
    # uses pmap; on a test with quadrature, this took about 3 minutes for the dataset
    # result = pmap(i-> likelihood(γ, β, ρ, α_0, α_1, α_2, t[i], x[i,:], z[i,:]), i for i = 1:N)
    
    # uses distributed for loop; on a test with quadrature, this took about 2 minutes for the dataset
    result = SharedArray{Float64}(N)

    if method == "quadrature"
        @sync @distributed for i = 1:N
            result[i] = likelihood_quadrature(γ, β, ρ, α_0, α_1, α_2, t[i], x[i,:], z[i,:])
        end

    elseif method == "ghk"
        @sync @distributed for i = 1:N
            result[i] = likelihood_ghk(γ, β, ρ, α_0, α_1, α_2, t[i], x[i,:], z[i,:])
        end

    elseif method == "accept_reject"
        @sync @distributed for i = 1:N
            result[i] = likelihood_accept_reject(γ, β, ρ, α_0, α_1, α_2, t[i], x[i,:], z[i,:])
        end
    end

    return result
end

function log_likelihood(γ::Array{Float64, 1}, β::Array{Float64, 1}, ρ::Float64, α_0::Float64, α_1::Float64, α_2::Float64,
    t::Array{Float64, 2}, x::Array{Float64, 2}, z::Array{Float64, 2})

    sum(log.(likelihood(γ, β, ρ, α_0, α_1, α_2, t, x, z)))

end

function log_likelihood(θ::Array{Float64, 1}, t::Array{Float64, 2}, x::Array{Float64, 2}, z::Array{Float64, 2})
    K_x = size(x)[2]
    K_z = size(z)[2]

    γ = θ[1:K_z]
    β = θ[(K_z+1):(K_x + K_z)]
    ρ = θ[K_x + K_z + 1]
    α_0 = θ[K_x + K_z + 2]
    α_1 = θ[K_x + K_z + 3]
    α_2 = θ[K_x + K_z + 4]

    log_likelihood(γ, β, ρ, α_0, α_1, α_2, t, x, z)
end 
